module solver
    implicit none

contains

subroutine compute_eig(Hk, eig, wavef)
    use interface, only: zheevr
    use constants, only: rank, prec, timer_debug, iband, nbands, nions
    implicit none
    complex(prec), intent(inout) :: Hk(:,:)
    real(prec), intent(out) :: eig(:)
    complex(prec), intent(out), optional :: wavef(:,:)
    
    character(len=1) :: jobz
    integer :: info, m
    complex(8), allocatable :: work8(:), z(:,:)
    real(8), allocatable :: rwork8(:)
    integer, allocatable :: iwork(:), isuppz(:)
    complex(8) :: workq(1)
    real(8) :: rworkq(1), vl, vu, abstol
    integer :: iworkq(1)
    integer :: lwork, lrwork, liwork
    real :: t_start, t_end

    if (timer_debug) call cpu_time(t_start)

    if (present(wavef)) then
        jobz = 'V'
        ! check the size of wavef 
        if (size(wavef, 1) /= nions .or. size(wavef, 2) /= nbands) then
            write(*,*) "[Error] wavef has wrong dimensions:", &
                       size(wavef, 1), "x", size(wavef, 2), &
                       ", expected:", nions, "x", nbands
            return
        end if
        allocate(z(nions, nbands))
    else
        jobz = 'N'
        allocate(z(1, 1))
    end if
    
    allocate(isuppz(2*nbands))

    abstol = 0.0_8
    vl = 0.0_8
    vu = 0.0_8
    
    ! workspace inquiry
    lwork = -1
    lrwork = -1
    liwork = -1
    
    call zheevr(jobz, 'I', 'U', nions, Hk, max(1, nions), vl, vu, iband(1), iband(2), abstol, &
                m, eig, z, max(1, nions), isuppz, workq, lwork, rworkq, lrwork, &
                iworkq, liwork, info)
    
    lwork = max(1, int(real(workq(1))))
    lrwork = max(1, int(rworkq(1)))
    liwork = max(1, iworkq(1))
    
    allocate(work8(lwork), rwork8(lrwork), iwork(liwork))
    
    ! real call
    call zheevr(jobz, 'I', 'U', nions, Hk, max(1, nions), vl, vu, iband(1), iband(2), abstol, &
                m, eig, z, max(1, nions), isuppz, work8, lwork, rwork8, lrwork, &
                iwork, liwork, info)
    
    if (timer_debug) then
        call cpu_time(t_end)
        write(*,"(12X,A,I4,1X,A,F8.3)") "[Eig] Rank", rank, "Call ZHEEVR finished!  Time cost (s): ", t_end-t_start
    end if
    
    if (info /= 0) then
        write(*,*) "[Error] Call ZHEEVR failed with info =", info
    end if
    
    ! output eigen functions
    if (present(wavef)) then
        wavef(1:nions, 1:nbands) = z(1:nions, 1:nbands)
    end if
    
    ! clean work arrays
    deallocate(work8, rwork8, iwork, isuppz, z)

end subroutine compute_eig

    subroutine calculate_k(kvec,eig,wavef)
        use constants, only: prec,nions
        use tbmodel, only: build_H
        real(prec), intent(in) :: kvec(3)
        real(prec), intent(out) :: eig(:)       
        complex(prec), intent(out), optional :: wavef(:,:)

        complex(prec), allocatable:: Hk(:,:)

        allocate(Hk(nions,nions))
        
        if (not(present(wavef))) then
            ! N for eigenvalues only , V also for eigenvectors
            call build_H(kvec,Hk)
            call compute_eig(Hk, eig)
        else
            call build_H(kvec,Hk)
            call compute_eig(Hk, eig, wavef)
        end if
        deallocate(Hk)    

    end subroutine calculate_k

    subroutine calculate_klist(klist,eig,wavef,tag)
        !  Fractional k-points are favored
        use constants, only: prec,nions,nbands
        use utils, only: k2frac_list
        real(prec), intent(in) :: klist(:,:)
        real(prec), intent(out), allocatable :: eig(:,:)
        complex(prec), intent(out), allocatable, optional :: wavef(:,:,:)
        character(len=1), intent(in), optional :: tag

        real(prec), allocatable :: kpoints(:,:)
        integer :: nkpts,i

        nkpts = size(klist,1)
        allocate(kpoints(nkpts,3))
        if (present(tag) .and. tag == 'C') then
            call k2frac_list(klist,kpoints)
        else 
            kpoints = klist
        end if

        if (.not. allocated(eig)) allocate(eig(nbands,nkpts))
        if (present(wavef) .and. .not. allocated(wavef)) allocate(wavef(nions,nbands,nkpts))

        do i = 1,nkpts
            if (present(wavef)) then
                call calculate_k(kpoints(i,:),eig(:,i),wavef(:,:,i))
            else
                call calculate_k(kpoints(i,:),eig(:,i))
            end if
        end do

    end subroutine calculate_klist


end module solver