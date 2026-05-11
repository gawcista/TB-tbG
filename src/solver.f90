module solver
    implicit none

contains

subroutine compute_eig(Hk, eig, wavef)
    use constants, only: rank, prec, timer_debug, iband, nbands, nions
    implicit none
    complex(prec), intent(inout) :: Hk(:,:)
    real(prec), intent(out) :: eig(:)
    complex(prec), intent(out), optional :: wavef(:,:)
    
    character(len=1) :: jobz
    complex(8), allocatable :: work8(:), z(:,:)
    real(8), allocatable :: rwork8(:)
    integer, allocatable :: iwork(:), isuppz(:)
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
    else
        jobz = 'N'
    end if

    call init_eig_workspace(jobz, work8, rwork8, iwork, isuppz, z)
    if (present(wavef)) then
        call compute_eig_cached(Hk, eig, jobz, work8, rwork8, iwork, isuppz, z, wavef)
    else
        call compute_eig_cached(Hk, eig, jobz, work8, rwork8, iwork, isuppz, z)
    end if
    
    if (timer_debug) then
        call cpu_time(t_end)
        write(*,"(12X,A,I4,1X,A,F8.3)") "[Eig] Rank", rank, "Call ZHEEVR finished!  Time cost (s): ", t_end-t_start
    end if

    deallocate(work8, rwork8, iwork, isuppz, z)

end subroutine compute_eig

    subroutine init_eig_workspace(jobz, work8, rwork8, iwork, isuppz, z)
        use interface, only: zheevr
        use constants, only: prec, iband, nbands, nions
        character(len=1), intent(in) :: jobz
        complex(8), allocatable, intent(out) :: work8(:), z(:,:)
        real(8), allocatable, intent(out) :: rwork8(:)
        integer, allocatable, intent(out) :: iwork(:), isuppz(:)

        complex(prec), allocatable :: Htmp(:,:)
        real(prec), allocatable :: eig_tmp(:)
        complex(8) :: workq(1)
        real(8) :: rworkq(1), vl, vu, abstol
        integer :: iworkq(1)
        integer :: info, m, lwork, lrwork, liwork

        allocate(Htmp(nions,nions), eig_tmp(nbands), isuppz(2*nbands))
        Htmp = (0.0_prec, 0.0_prec)
        if (jobz == 'V') then
            allocate(z(nions, nbands))
        else
            allocate(z(1, 1))
        end if

        abstol = 0.0_8
        vl = 0.0_8
        vu = 0.0_8
        lwork = -1
        lrwork = -1
        liwork = -1

        call zheevr(jobz, 'I', 'U', nions, Htmp, max(1, nions), vl, vu, iband(1), iband(2), abstol, &
                    m, eig_tmp, z, size(z,1), isuppz, workq, lwork, rworkq, lrwork, &
                    iworkq, liwork, info)

        lwork = max(1, int(real(workq(1))))
        lrwork = max(1, int(rworkq(1)))
        liwork = max(1, iworkq(1))
        allocate(work8(lwork), rwork8(lrwork), iwork(liwork))

        deallocate(Htmp, eig_tmp)
    end subroutine init_eig_workspace

    subroutine compute_eig_cached(Hk, eig, jobz, work8, rwork8, iwork, isuppz, z, wavef)
        use interface, only: zheevr
        use constants, only: prec, iband, nbands, nions
        complex(prec), intent(inout) :: Hk(:,:)
        real(prec), intent(out) :: eig(:)
        character(len=1), intent(in) :: jobz
        complex(8), intent(inout) :: work8(:), z(:,:)
        real(8), intent(inout) :: rwork8(:)
        integer, intent(inout) :: iwork(:), isuppz(:)
        complex(prec), intent(out), optional :: wavef(:,:)

        real(8) :: vl, vu, abstol
        integer :: info, m

        abstol = 0.0_8
        vl = 0.0_8
        vu = 0.0_8

        call zheevr(jobz, 'I', 'U', nions, Hk, max(1, nions), vl, vu, iband(1), iband(2), abstol, &
                    m, eig, z, size(z,1), isuppz, work8, size(work8), rwork8, size(rwork8), &
                    iwork, size(iwork), info)

        if (info /= 0) then
            write(*,*) "[Error] Call ZHEEVR failed with info =", info
        end if

        if (present(wavef)) then
            wavef(1:nions, 1:nbands) = z(1:nions, 1:nbands)
        end if
    end subroutine compute_eig_cached

    subroutine calculate_k(kvec,eig,wavef)
        use constants, only: prec,nions
        use tbmodel, only: build_H
        real(prec), intent(in) :: kvec(3)
        real(prec), intent(out) :: eig(:)       
        complex(prec), intent(out), optional :: wavef(:,:)

        complex(prec), allocatable:: Hk(:,:)

        allocate(Hk(nions,nions))
        
        if (.not. present(wavef)) then
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
        use tbmodel, only: build_H
        real(prec), intent(in) :: klist(:,:)
        real(prec), intent(out), allocatable :: eig(:,:)
        complex(prec), intent(out), allocatable, optional :: wavef(:,:,:)
        character(len=1), intent(in), optional :: tag

        real(prec), allocatable :: kpoints(:,:)
        complex(prec), allocatable :: Hk(:,:)
        complex(8), allocatable :: work8(:), z(:,:)
        real(8), allocatable :: rwork8(:)
        integer, allocatable :: iwork(:), isuppz(:)
        character(len=1) :: jobz
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

        if (present(wavef)) then
            jobz = 'V'
        else
            jobz = 'N'
        end if
        call init_eig_workspace(jobz, work8, rwork8, iwork, isuppz, z)
        allocate(Hk(nions,nions))

        do i = 1,nkpts
            call build_H(kpoints(i,:),Hk)
            if (present(wavef)) then
                call compute_eig_cached(Hk,eig(:,i),jobz,work8,rwork8,iwork,isuppz,z,wavef(:,:,i))
            else
                call compute_eig_cached(Hk,eig(:,i),jobz,work8,rwork8,iwork,isuppz,z)
            end if
        end do

        deallocate(Hk, work8, rwork8, iwork, isuppz, z)

    end subroutine calculate_klist


end module solver
