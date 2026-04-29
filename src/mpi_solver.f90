module mpi_solver
    use mpi_f08
    use tbmodel


    implicit none
    
contains

    subroutine init_mpi()
        use constants, only: rank,nproc,init_constants
        use parser, only: input_parser
        integer :: ierr
        character(len=64) :: input_file

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
        call init_constants()
        if (command_argument_count() >= 1) then
            call get_command_argument(1, input_file)
            call input_parser(input_file)
        else
            call input_parser()
        end if
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
        integer :: base_load, remainder, i, cache_size

        if (.not. allocated(counts)) allocate(counts(nprocs))
        if (.not. allocated(displs)) allocate(displs(nprocs))
        cache_size = get_cache_size()
        
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

        if (present(ngroups)) then
            if (.not. allocated(ngroups)) allocate(ngroups(nprocs))
            do i = 1, nprocs
                ngroups(i) = (counts(i) + cache_size - 1) / cache_size
            end do
        end if


    end subroutine calculate_workload

    integer function get_cache_size() result(cache_size)
        use constants, only: ncache

        if (ncache > 0) then
            cache_size = ncache
        else
            cache_size = 1
        end if
    end function get_cache_size

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
        integer :: ngroups, maxgroups, igroup, group_start, group_end, group_size, cache_size
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
        cache_size = get_cache_size()
        
        do igroup = 1, maxgroups
            if (timer_mpi) call cpu_time(t_mpi)
            if (igroup <= ngroups) then
                group_start = (igroup - 1) * cache_size + 1
                group_end = min(igroup * cache_size, nkpts_local)
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
                        i_group_start_local = (igroup - 1) * cache_size + 1
                        i_group_end_local = min(igroup * cache_size, i_nkpts_local)
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
    use ioutils, only: readKPOINTS,writeEELS,writeEELSmatrix
    use constants, only: f_kpoint,prec,basis_rec,rank,nproc,ncache,nomega,omegalist,timer_cpu,timer_mpi,nbands,eels_mode,iband,write_matrix,f_matrix
    use utils,only: cross2d,k2cart_list,k2cart_list,k2frac_list
    use eels, only: calculate_IPF_klist,calculate_eels
    real(prec),intent(in) :: q_list(:,:)
    real(prec), allocatable :: eloss(:,:)
    character(1), intent(in), optional :: tag
    real(prec),allocatable :: klist(:,:),klist_local(:,:),klist_group(:,:),q_list_frac(:,:),q_list_cart(:,:)
    integer,allocatable :: counts(:), displs(:), ngroups_all(:)
    integer :: nkpts,nqpts,nkpts_local,ik_start,ik_end,ngroups,ierr
    integer :: igroup,group_start,group_end,group_size
    integer :: cache_size
    complex(prec),allocatable :: IPF(:,:),IPF_local(:,:),IPF_group(:,:)
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
    allocate(IPF_local(nomega,nqpts))
    IPF_local = (0.0_prec, 0.0_prec)
    if (write_matrix) then
        allocate(matrix_contrib_local(nbands, nbands, nomega, nqpts))
        matrix_contrib_local = (0.0_prec, 0.0_prec)
    end if

    ngroups = ngroups_all(rank + 1)
    cache_size = get_cache_size()
    write(*, '(4X,A,I4,1X,A,I4,1X,A,1X,I4,1X,A)') "[EELS] Rank", rank,"calculating", nkpts_local, "kpoints in", ngroups, "groups"

    do igroup = 1, maxval(ngroups_all)
        if (timer_mpi) call cpu_time(t_mpi)
        ! working part
        if (igroup <= ngroups) then
            ! get group kpoints from local 
            group_start = (igroup - 1) * cache_size + 1
            group_end = min(igroup * cache_size, nkpts_local)
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
            IPF_local = IPF_local + IPF_group
            deallocate(klist_group)
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

    call MPI_Reduce(IPF_local, IPF, nomega*nqpts, MPI_DOUBLE_COMPLEX, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

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
            call writeEELSmatrix(q_list_cart, IPF, matrix_contrib_global, dS)
        end if
    end if

    deallocate(IPF_local)

end subroutine calculate_eels_mpi

    subroutine calculate_dos_mpi()
        use constants, only: prec,rank,kpt_select,f_kpoint,nk_path,nedos,erange,sigma_dos,gaussian_cutoff,f_dos
        use ioutils, only: readKPOINTS,writeDOS
        use utils, only: generate_mesh_gamma, generate_mesh_mp
        real(prec), allocatable :: klist_frac(:,:), dos(:)
        real(prec), allocatable :: energy_grid(:)
        integer :: ierr,nkpts,i,nedos_local

        if (rank == 0) then
            write(*,"(A)") "[DOS] Starting DOS calculation"
            if (nedos < 2) then
                write(*,"(A)") "[DOS] Error: NEDOS must be greater than 1"
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            if (sigma_dos <= 0.0_prec) then
                write(*,"(A)") "[DOS] Error: SIGMA must be positive"
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            if (erange(2) <= erange(1)) then
                write(*,"(A)") "[DOS] Error: ERANGE upper bound must be greater than lower bound"
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            ! Generate energy grid
            allocate(energy_grid(nedos))
            do i = 1, nedos
                energy_grid(i) = erange(1) + (erange(2)-erange(1))*(i-1)/(nedos-1)
            end do

            ! Read k-point grid
            if (kpt_select) then
                call readKPOINTS(f_kpoint, klist_frac)
            else
                ! Default: use Gamma-centered uniform mesh
                ! If KPOINTS file doesn't exist, generate default mesh
                ! Try reading file first, generate default mesh if fails
                call readKPOINTS(f_kpoint, klist_frac)
            end if
            if (.not. allocated(klist_frac)) then
                write(*,"(A)") "[DOS] Error: no k-points were read"
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            nkpts = size(klist_frac,1)
            if (nkpts <= 0) then
                write(*,"(A)") "[DOS] Error: no k-points were read"
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            write(*,"(A,I8)") "[DOS] Number of k-points: ", nkpts
        end if

        ! Broadcast energy grid and k-point info
        nedos_local = nedos
        if (rank == 0) write(*,"(A,I8)") "[DOS] Broadcasting nedos = ", nedos_local
        call MPI_Bcast(nedos_local, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (rank /= 0) allocate(energy_grid(nedos_local))
        call MPI_Bcast(energy_grid, nedos_local, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(nkpts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (.not. allocated(klist_frac)) allocate(klist_frac(nkpts,3))
        call MPI_Bcast(klist_frac, nkpts*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! Call core DOS calculation function
        call calculate_dos_klist_mpi(klist_frac, energy_grid, dos)

        if (rank == 0) then
            ! Output DOS results
            call writeDOS(energy_grid, dos, f_dos)
            write(*,"(A)") "[DOS] DOS calculation completed"
        end if

        deallocate(energy_grid)
        if (allocated(klist_frac)) deallocate(klist_frac)
        if (allocated(dos)) deallocate(dos)
    end subroutine calculate_dos_mpi

    subroutine calculate_dos_klist_mpi(klist, energy_grid, dos)
        use constants, only: prec,ncache,nbands,nions,nproc,rank,timer_mpi,timer_cpu,sigma_dos,gaussian_cutoff,dos_normalize
        use utils, only: k2frac_list, add_gaussian_to_dos
        use solver, only: calculate_klist
        integer :: ierr
        real(prec), intent(in) :: klist(:,:)
        real(prec), intent(in) :: energy_grid(:)
        real(prec), intent(out), allocatable :: dos(:)

        integer :: nkpts,nkpts_local,ik_start,ik_end
        real :: t_start,t_end, t_mpi
        integer,allocatable :: counts(:), displs(:), ngroups_all(:)
        real(prec), allocatable :: eig_group(:,:)
        real(prec), allocatable :: klist_local(:,:), klist_group(:,:)
        integer :: ngroups, maxgroups, igroup, group_start, group_end, group_size, cache_size
        integer :: i, j, iband
        type(MPI_Status) :: mpi_status
        integer :: i_nkpts_local, i_ik_start, i_group_start_local, i_group_end_local, i_group_size, i_group_start_global
        real(prec), allocatable :: dos_local(:), dos_recv(:)
        real(prec), allocatable :: eig_recv(:,:)

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

        ! Initialize local DOS array
        allocate(dos_local(size(energy_grid)))
        dos_local = 0.0_prec

        ! Convert k-point coordinates (if needed)
        allocate(klist_local(nkpts_local,3))
        klist_local = klist(ik_start:ik_end,:)

        write(*, '(4X,A,1X,I4,1X,A,I4)') "[DOS]", nkpts_local, "k-points was assigned to rank", rank

        if (rank == 0) then
            allocate(dos(size(energy_grid)))
            dos = 0.0_prec
        end if

        call cpu_time(t_start)

        ngroups = ngroups_all(rank+1)
        maxgroups = maxval(ngroups_all)
        cache_size = get_cache_size()

        do igroup = 1, maxgroups
            if (timer_mpi) call cpu_time(t_mpi)
            if (igroup <= ngroups) then
                group_start = (igroup - 1) * cache_size + 1
                group_end = min(igroup * cache_size, nkpts_local)
                group_size = group_end - group_start + 1

                allocate(klist_group(group_size, 3))
                klist_group = klist_local(group_start:group_end, :)

                allocate(eig_group(nbands, group_size))
                call calculate_klist(klist_group, eig_group)

                ! Add Gaussian contribution from each eigenvalue to local DOS
                do i = 1, group_size
                    do iband = 1, nbands
                        call add_gaussian_to_dos(eig_group(iband, i), sigma_dos, gaussian_cutoff, energy_grid, dos_local)
                    end do
                end do

                deallocate(eig_group, klist_group)
            end if

            if (timer_mpi .and. igroup <= ngroups) then
                call cpu_time(t_end)
                write(*, '(4X,A,1X,I4,1X,A,I3,A,I3,1X,A,F12.3)') &
                    "[DOS] Rank", rank, "group", igroup, "/", ngroups, &
                    "finished! CPU time:", t_end-t_mpi
            end if
        end do

        ! Use MPI_Reduce to gather DOS contributions from all processes
        call MPI_Reduce(dos_local, dos, size(energy_grid), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        ! Apply k-point normalization on rank 0 if requested
        if (rank == 0 .and. dos_normalize) then
            dos = dos / real(nkpts, prec)
            write(*,"(A,I8)") "[DOS] Applied k-point normalization, Nk = ", nkpts
        end if

        if (timer_cpu)  then
            call cpu_time(t_end)
            write(*, '(4X,A,1X,I4,1X,A,F12.3)') "[DOS] Calculationg on rank",rank,&
                                             'finished! CPU time:', t_end-t_start
        end if

        deallocate(klist_local, counts, displs, ngroups_all, dos_local)
        if (rank /= 0 .and. allocated(dos)) deallocate(dos)
    end subroutine calculate_dos_klist_mpi

end module mpi_solver
