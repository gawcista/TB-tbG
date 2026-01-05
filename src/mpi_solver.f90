module mpi_solver
    use mpi_f08
    use tbmodel


    implicit none
    
contains

    subroutine init_mpi()
        use constants, only: rank,nproc,init_constants
        use parser, only: input_parser
        integer :: ierr

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
        call init_constants()
        call input_parser()
    end subroutine init_mpi

    subroutine init_POSCAR_mpi(filename)
        use constants, only: rank,nions,nbands,iband,basis,basis_rec,position_frac,position_cart,f_poscar,timer_mpi
        use ioutils, only: parsePOSCAR
        integer :: ierr
        character(len=*), intent(in), optional :: filename
        real :: t_start,t_end


        if (rank == 0) then
            if (present(filename)) then
                call parsePOSCAR(filename)
            else
                call parsePOSCAR(f_poscar)
            end if
        end if

        ! Broadcast the lattice structure to all processes
        if (rank == 0 .and. timer_mpi) then
            call cpu_time(t_start)
            write(*, '(4X,A,1X,$)') "[MPI] Updating POSCAR data to all MPI processes ..."
        end if

        call MPI_Bcast(nions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(nbands, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(iband, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(basis, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(basis_rec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (.not. allocated(position_frac)) then
            allocate(position_frac(nions,3))
        end if
        if (.not. allocated(position_cart)) then
            allocate(position_cart(nions,3))  
        end if
        call MPI_Bcast(position_frac, size(position_frac), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(position_cart, size(position_cart), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (rank == 0 .and. timer_mpi) then
            call cpu_time(t_end)
            write(*, '(A,1X,F12.3)') "Done! Time elapsed (s): ", t_end - t_start
        end if

    end subroutine init_POSCAR_mpi

    subroutine calculate_workload(total_size, nprocs, counts, displs, ngroups)
        ! Assign MPI jobs
        use constants, only: ncache
        integer, intent(in) :: total_size, nprocs
        integer, intent(out), allocatable :: counts(:), displs(:)
        integer, intent(out), allocatable, optional :: ngroups(:)
        integer :: base_load, remainder, i, n_groups

        if (.not. allocated(counts)) allocate(counts(nprocs))
        if (.not. allocated(displs)) allocate(displs(nprocs))
        
        base_load = total_size / nprocs
        remainder = mod(total_size, nprocs)
        
        ! Initialize all with base load
        counts = base_load
        
        ! Distribute remainder tasks starting from rank 1 (index 2) to avoid giving extra to rank 0
        if (remainder > 0) then
            ! Rank 0 (index 1) gets base_load, rank 1 to rank remainder (indices 2 to remainder+1) get base_load + 1
            ! Ensure we don't go beyond array bounds
            do i = 2, remainder+1
                counts(i) = base_load + 1
            end do
        end if
        
        displs(1) = 0
        do i = 2, nprocs
            displs(i) = displs(i-1) + counts(i-1)
        end do

        if (present(ngroups) .and. ncache>0) then
            if (.not. allocated(ngroups)) allocate(ngroups(nprocs))
            do i = 1, nprocs
                n_groups = counts(i) / ncache
                if (mod(counts(i), ncache) /= 0) then
                    ngroups(i) = n_groups + 1
                else
                    ngroups(i) = n_groups
                end if
            end do
        end if


    end subroutine calculate_workload

subroutine calculate_klist_mpi(klist,eig,wavef,tag)
        use constants, only: prec,ncache,nbands,nions,nproc,rank,timer_mpi,timer_cpu
        use utils, only: k2frac_list
        use solver, only: calculate_klist
        integer :: ierr
        real(prec), intent(in) :: klist(:,:)
        real(prec), intent(out), allocatable :: eig(:,:)
        complex(prec), intent(out), allocatable, optional :: wavef(:,:,:)
        character(len=1), intent(in), optional :: tag

        integer :: nkpts,nkpts_local,ik_start,ik_end
        real :: t_start,t_end, t_mpi
        integer,allocatable :: counts(:), displs(:), ngroups_all(:)
        real(prec), allocatable :: eig_group(:,:)
        complex(prec), allocatable :: wavef_group(:,:,:)
        real(prec), allocatable :: klist_local(:,:), klist_group(:,:)
        integer :: ngroups, maxgroups, igroup, group_start, group_end, group_size
        integer :: i, j
        type(MPI_Status) :: mpi_status
        integer :: i_nkpts_local, i_ik_start, i_group_start_local, i_group_end_local, i_group_size, i_group_start_global
        real(prec), allocatable :: eig_recv(:,:)
        complex(prec), allocatable :: wavef_recv(:,:,:)
        
        nkpts = size(klist,1)
        allocate(counts(nproc), displs(nproc), ngroups_all(nproc))
        if (rank == 0) then
            call calculate_workload(nkpts, nproc, counts, displs,ngroups_all)
        end if
        call MPI_Bcast(counts, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(displs, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ngroups_all, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
        
        nkpts_local = counts(rank+1)
        ik_start = displs(rank+1) + 1
        ik_end = ik_start + nkpts_local - 1

        if (.not. present(tag) .or. tag == 'D') then
            klist_local = klist(ik_start:ik_end,:)
        else
            allocate(klist_local(nkpts_local,3))
            call k2frac_list(klist(ik_start:ik_end,:),klist_local)
        end if
        write(*, '(4X,A,1X,I4,1X,A,I4)') "[Solver]", nkpts_local, "k-points was assigned to rank", rank

        if (rank == 0) then
            allocate(eig(nbands, nkpts))
            if (present(wavef)) then
                allocate(wavef(nions, nbands, nkpts))
            end if
        end if

        call cpu_time(t_start)

        ngroups = ngroups_all(rank+1)
        maxgroups = maxval(ngroups_all)
        
        do igroup = 1, maxgroups
            if (timer_mpi) call cpu_time(t_mpi)
            if (igroup <= ngroups) then
                group_start = (igroup - 1) * ncache + 1
                group_end = min(igroup * ncache, nkpts_local)
                group_size = group_end - group_start + 1
                
                allocate(klist_group(group_size, 3))
                klist_group = klist_local(group_start:group_end, :)
                
                if (present(wavef)) then
                    allocate(wavef_group(nions, nbands, group_size))
                    allocate(eig_group(nbands, group_size))
                    call calculate_klist(klist_group, eig_group, wavef_group)
                else
                    allocate(eig_group(nbands, group_size))
                    call calculate_klist(klist_group, eig_group)
                end if
            end if
            
            if (rank /= 0) then
                if (igroup <= ngroups) then
                    call MPI_Send(eig_group, nbands*group_size, MPI_DOUBLE_PRECISION, &
                                  0, igroup, MPI_COMM_WORLD, ierr)
                    if (present(wavef)) then
                        call MPI_Send(wavef_group, nions*nbands*group_size, &
                                      MPI_DOUBLE_COMPLEX, 0, igroup+nproc, MPI_COMM_WORLD, ierr)
                    end if
                end if
            else
                do i = 1, nproc-1
                    if (igroup <= ngroups_all(i+1)) then
                        i_nkpts_local = counts(i+1)
                        i_ik_start = displs(i+1) + 1
                        i_group_start_local = (igroup - 1) * ncache + 1
                        i_group_end_local = min(igroup * ncache, i_nkpts_local)
                        i_group_size = i_group_end_local - i_group_start_local + 1
                        i_group_start_global = i_ik_start + i_group_start_local - 1
                        
                        allocate(eig_recv(nbands, i_group_size))
                        call MPI_Recv(eig_recv, nbands*i_group_size, MPI_DOUBLE_PRECISION, &
                                      i, igroup, MPI_COMM_WORLD, mpi_status, ierr)
                        
                        do j = 1, i_group_size
                            eig(:, i_group_start_global + j - 1) = eig_recv(:, j)
                        end do
                        deallocate(eig_recv)
                        
                        if (present(wavef)) then
                            allocate(wavef_recv(nions, nbands, i_group_size))
                            call MPI_Recv(wavef_recv, nions*nbands*i_group_size, &
                                          MPI_DOUBLE_COMPLEX, i, igroup+nproc, &
                                          MPI_COMM_WORLD, mpi_status, ierr)
                            
                            do j = 1, i_group_size
                                wavef(:, :, i_group_start_global + j - 1) = wavef_recv(:, :, j)
                            end do
                            deallocate(wavef_recv)
                        end if
                    end if
                end do
                
                if (igroup <= ngroups) then
                    do j = 1, group_size
                        eig(:, ik_start + group_start + j - 2) = eig_group(:, j)
                        if (present(wavef)) then
                            wavef(:, :, ik_start + group_start + j - 2) = wavef_group(:, :, j)
                        end if
                    end do
                end if
            end if
            
            if (igroup <= ngroups) then
                deallocate(klist_group, eig_group)
                if (allocated(wavef_group)) deallocate(wavef_group)
            end if
            
            if (timer_mpi .and. igroup <= ngroups) then
                call cpu_time(t_end)
                write(*, '(4X,A,1X,I4,1X,A,I3,A,I3,1X,A,F12.3)') &
                    "[Solver] Rank", rank, "group", igroup, "/", ngroups, &
                    "finished! CPU time:", t_end-t_mpi
            end if
        end do

        if (timer_cpu)  then
            call cpu_time(t_end)
            write(*, '(4X,A,1X,I4,1X,A,F12.3)') "[Solver] Calculationg on rank",rank,&
                                             'finished! CPU time:', t_end-t_start
        end if

        deallocate(klist_local, counts, displs, ngroups_all)

end subroutine calculate_klist_mpi

subroutine calculate_band_mpi(kpath, nk)
    use constants, only: prec,kpath_default,rank,kpt_select,f_kpoint,nk_path
    use ioutils, only: readKPOINTS,writeBand,writeLabels
    real(prec), allocatable, optional:: kpath(:,:,:)
    integer, optional :: nk
    real(prec), allocatable:: paths(:,:,:),klist_frac(:,:),klist_cart(:,:),xlist(:),lable_positions(:),eig(:,:)
    real(prec), allocatable:: klist_tmp(:,:)
    character(len=3), allocatable :: labels(:)
    integer :: ierr,nkpts,nk_local,i,npath
    
    if (rank == 0) then
        if (present(kpath)) then
            paths = kpath
        else
            if (kpt_select) then
                call readKPOINTS(f_kpoint,klist_tmp)
                npath = size(klist_tmp,1)/2
                allocate(paths(npath,2,3))
                do i = 1, npath
                    paths(i,1,:) = klist_tmp(2*i-1,:)
                    paths(i,2,:) = klist_tmp(2*i,:)
                end do
            else
                call kpath_default(paths)
            end if
        end if
        if (present(nk)) then 
            nk_local = nk
        else 
            nk_local = nk_path
        end if
        call generate_kpath(paths, nk_local, klist_frac, klist_cart, xlist, labels, lable_positions)
        nkpts = size(klist_frac,1)
    end if


    call MPI_Bcast(nkpts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if ( .not. allocated(klist_frac)) allocate(klist_frac(nkpts,3))
    call MPI_Bcast(klist_frac, nkpts*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call calculate_klist_mpi(klist_frac,eig)

    if (rank == 0) then
        call writeBand(xlist, eig)
        call writeLabels(labels, lable_positions)
    end if
end subroutine calculate_band_mpi

subroutine calculate_eels_mpi(q_list,tag)
    ! q by default: Cartesian
    use ioutils, only: readKPOINTS,writeEELS
    use constants, only: f_kpoint,prec,basis_rec,rank,nproc,ncache,nomega,omegalist,timer_cpu,timer_mpi,nbands,eels_mode,iband,write_matrix,f_matrix
    use utils,only: cross2d,k2cart_list,k2cart_list,k2frac_list
    use eels, only: calculate_IPF_klist,calculate_eels
    type(MPI_Status) :: mpi_status
    real(prec),intent(in) :: q_list(:,:)
    real(prec), allocatable :: eloss(:,:)
    character(1), intent(in), optional :: tag
    real(prec),allocatable :: klist(:,:),klist_local(:,:),klist_group(:,:),q_list_frac(:,:),q_list_cart(:,:)
    integer,allocatable :: counts(:), displs(:), ngroups_all(:)
    integer :: i,nkpts,nqpts,nkpts_local,ik_start,ik_end,ngroups,ierr
    integer :: igroup,group_start,group_end,group_size
    complex(prec),allocatable :: IPF(:,:),IPF_group(:,:),IPF_recv(:,:)
    real(prec) :: dS,t_start,t_end,t_mpi
    complex(prec), allocatable :: matrix_contrib_local(:,:,:,:), matrix_contrib_global(:,:,:,:)
    complex(prec), allocatable :: matrix_contrib_group(:,:,:,:)

    nqpts = size(q_list,1)
    allocate(q_list_frac(nqpts,3),q_list_cart(nqpts,3))
    allocate(counts(nproc),displs(nproc),ngroups_all(nproc))
    if (rank == 0) then
        if (timer_cpu) call cpu_time(t_start)
        if (present(tag) .and. tag=='D') then
            q_list_frac = q_list
            call k2cart_list(q_list,q_list_cart)
        else
            call k2frac_list(q_list,q_list_frac)
            q_list_cart = q_list
        end if
        call readKPOINTS(f_kpoint,klist)
        nkpts = size(klist,1)
        dS = cross2d(basis_rec(1,:),basis_rec(2,:))/nkpts
        call calculate_workload(nkpts,nproc,counts,displs,ngroups_all)
        write(*,"(A,I4,1X,A,I4,1X,A)") "[Main] Calculating EELS of ",nqpts,"q-points on ",nkpts,"k-points"
        select case (eels_mode)
            case(0)
                write(*,'(4X,A)') '[EELS] EELS calulation mode is set to 0:ALL'
            case(1)
                write(*,'(4X,A)') '[EELS] EELS calulation mode is set to 1:intraband only'
            case(2)
                write(*,'(4X,A)') '[EELS] EELS calulation mode is set to 2:interband only'
        end select
        write(*,"(4X,A,I6,1X,A,I4,1X,A,I4,1X,A)") "[EELS] ",nbands," bands (#",iband(1),"-",iband(2),") will be considered in the calculation."
        write(*, '(4X,A,I4,1X,A,F5.3,A,F5.3,A)') "[EELS]",nomega,"omega points are set in [",omegalist(1),",",omegalist(nomega),"]"
        ! Output matrix elements if requested
        if (write_matrix) then
            write(*,'(4X,A)') '[EELS] Writing transition matrix contributions to file: '//trim(f_matrix)
        end if
    end if
    ! Broadcast nkpts and workload
    call MPI_Bcast(nkpts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(q_list_frac, nqpts * 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(q_list_cart, nqpts * 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(counts, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(displs, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ngroups_all, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    nkpts_local = counts(rank + 1)
    ik_start = displs(rank + 1) + 1
    ik_end = ik_start + nkpts_local - 1
    if (rank/=0) allocate(klist(nkpts,3))
    call MPI_Bcast(klist, nkpts * 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    ! obtain local jobs
    allocate(klist_local(nkpts_local, 3))
    klist_local = klist(ik_start:ik_end, :)
    if (rank /= 0) deallocate(klist)
    
    if (rank ==0) then
        allocate(IPF(nomega,nqpts))
        IPF = (0.0_prec, 0.0_prec)
        if (write_matrix) then
            allocate(matrix_contrib_global(nbands, nbands, nomega, nqpts))
            matrix_contrib_global = (0.0_prec, 0.0_prec)
        end if
    end if
    if (write_matrix) then
        allocate(matrix_contrib_local(nbands, nbands, nomega, nqpts))
        matrix_contrib_local = (0.0_prec, 0.0_prec)
    end if

    ngroups = ngroups_all(rank + 1)
    write(*, '(4X,A,I4,1X,A,I4,1X,A,1X,I4,1X,A)') "[EELS] Rank", rank,"calculating", nkpts_local, "kpoints in", ngroups, "groups"

    do igroup = 1, maxval(ngroups_all)
        if (timer_mpi) call cpu_time(t_mpi)
        ! working part
        if (igroup <= ngroups) then
            ! get group kpoints from local 
            group_start = (igroup - 1) * ncache + 1
            group_end = min(igroup * ncache, nkpts_local)
            group_size = group_end - group_start + 1
            allocate(klist_group(group_size, 3))
            klist_group = klist_local(group_start:group_end, :)
            ! calculate IPF
            allocate(IPF_group(nomega,nqpts))
            IPF_group = (0.0_prec, 0.0_prec)
            if (write_matrix) then
                allocate(matrix_contrib_group(nbands, nbands, nomega, nqpts))
                matrix_contrib_group = (0.0_prec, 0.0_prec)
                call calculate_IPF_klist(omegalist, klist_group, q_list_frac, IPF_group, matrix_contrib_group)
                ! accumulate local matrix contributions
                matrix_contrib_local = matrix_contrib_local + matrix_contrib_group
                deallocate(matrix_contrib_group)
            else
                call calculate_IPF_klist(omegalist, klist_group, q_list_frac, IPF_group)
            end if
            deallocate(klist_group)
        end if
        
        if (rank /= 0) then
            ! send result to rank0
            if (igroup <= ngroups) call MPI_Send(IPF_group, nomega * nqpts, MPI_DOUBLE_COMPLEX, 0, igroup, MPI_COMM_WORLD, ierr)
        else ! rank0: recieve
            do i = 1, nproc - 1
                if (igroup <= ngroups_all(i + 1)) then
                    allocate(IPF_recv(nomega,nqpts))
                    call MPI_Recv(IPF_recv, nomega*nqpts, MPI_DOUBLE_COMPLEX, i, igroup, MPI_COMM_WORLD, mpi_status, ierr)
                    IPF = IPF + IPF_recv
                    deallocate(IPF_recv)
                end if
            end do
            ! add local IPF calculated on rank0
            if (igroup <= ngroups) IPF = IPF + IPF_group
        end if

        if (allocated(IPF_group)) deallocate(IPF_group)

        if (timer_mpi .and. igroup <= ngroups) then
            call cpu_time(t_end)
            write(*, '(4X,A,1X,I4,1X,A,I3,A,I3,1X,A,F12.3)') &
                "[EELS] Rank", rank, "group", igroup, "/", ngroups, &
                "finished! Time:", t_end - t_mpi
        end if
    end do
    
    deallocate(klist_local, counts, displs, ngroups_all)

    ! Reduce matrix contributions from all processes if needed
    if (write_matrix) then
        call MPI_Reduce(matrix_contrib_local, matrix_contrib_global, &
                        nbands*nbands*nomega*nqpts, MPI_DOUBLE_COMPLEX, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

    if (rank == 0) then
        call calculate_eels(IPF*dS,q_list_cart,eloss)
        write(*,"(A,$)") "[Main] EELS calculation finished!"
        if (timer_cpu) then
            call cpu_time(t_end)
            write(*,"(1X,A,F12.3)") 'CPU time:', t_end-t_start
        else 
            write(*,*)
        end if
        call writeEELS(q_list_cart,omegalist,eloss)
        ! Write matrix contributions if requested
        if (write_matrix) then
            call write_matrix_contributions(q_list_cart, omegalist, IPF, matrix_contrib_global, dS)
        end if
    end if

end subroutine calculate_eels_mpi

subroutine write_matrix_contributions(q_list_cart, omegalist, IPF_total, matrix_contrib, dS)
    use constants, only: prec, nbands, f_matrix, iband
    use eels, only: V_q
    real(prec), intent(in) :: q_list_cart(:,:), omegalist(:)
    complex(prec), intent(in) :: IPF_total(:,:)  ! (nomega, nqpts)
    complex(prec), intent(in) :: matrix_contrib(:,:,:,:)  ! (nbands, nbands, nomega, nqpts)
    real(prec), intent(in) :: dS
    integer :: iq, iomega, i, j, unit, nqpts, nomega
    complex(prec) :: Vq, f_pi_ij, epsilon_ij, IPF_ij
    complex(prec) :: epsilon_total, IPF_total_value
    real(prec) :: q_vec(3), omega_val
    integer :: ib_start, ib_end
    
    ib_start = iband(1)
    ib_end = iband(2)
    
    nqpts = size(q_list_cart, 1)
    nomega = size(omegalist)
    
    unit = 500
    open(unit, file=f_matrix, status='replace', action='write')
    
    do iq = 1, nqpts
        q_vec = q_list_cart(iq, :)
        Vq = V_q(q_vec)
        do iomega = 1, nomega
            omega_val = omegalist(iomega)
            ! Calculate total IPF for this (q, omega)
            epsilon_total = 1.0_prec - Vq * IPF_total(iomega, iq)
            IPF_total_value = -imag(1.0_prec / epsilon_total)
            
            ! Write comment line
            write(unit, '(A,3F12.8,A,F12.8,A,F12.8)') '# q = [', q_vec, '] omega = ', omega_val, ' IPF_total = ', real(IPF_total_value)
            
            ! Write each transition contribution
            do j = ib_start, ib_end   ! iband_k
                do i = ib_start, ib_end   ! iband_kq
                    ! Note: matrix indices are (iband_kq, iband_k)
                    f_pi_ij = matrix_contrib(i, j, iomega, iq) * dS
                    epsilon_ij = 1.0_prec - Vq * f_pi_ij
                    IPF_ij = -imag(1.0_prec / epsilon_ij)
                    write(unit, '(2I6,1X,F20.12)') i, j, real(IPF_ij)
                end do
            end do
            
            ! Blank line between blocks
            write(unit, *)
        end do
    end do
    
    close(unit)
    write(*, '(4X,A)') '[EELS] Transition matrix contributions written to file: '//trim(f_matrix)
    
end subroutine write_matrix_contributions

end module mpi_solver
