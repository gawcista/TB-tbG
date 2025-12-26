module solver
    implicit none

contains

    subroutine compute_eig_full(Hk, eig, jobz)
        use interface, only: zheevd
        use constants, only: rank, prec, timer_debug
        complex(prec), intent(inout) :: Hk(:,:)
        real(prec), intent(out) :: eig(:)
        character(len=1), intent(in) :: jobz

        integer :: n, lda, info
        complex(8), allocatable :: work8(:)
        real(8), allocatable :: rwork8(:)
        integer, allocatable :: iwork(:)
        complex(8) :: workq(1)
        real(8) :: rworkq(1)
        integer :: iworkq(1)
        integer :: lwork, lrwork, liwork
        real :: t_start,t_end

        if (timer_debug) call cpu_time(t_start)

        ! jobz is required by caller; pass it directly to zheevd
        n = size(Hk,1)
        lda = max(1,n)

        lwork = -1
        lrwork = -1
        liwork = -1
        ! trial call
        call zheevd(jobz,'U', n, Hk, lda, eig, workq, lwork, rworkq, lrwork, iworkq, liwork, info)
        lwork = max(1, int(real(workq(1))))
        lrwork = max(1, int(rworkq(1)))
        liwork = max(1, iworkq(1))

        allocate(work8(lwork),rwork8(lrwork),iwork(liwork))

        call zheevd(jobz,'U', n, Hk, lda, eig, work8, lwork, rwork8, lrwork, iwork, liwork, info)

        if (timer_debug) then
            call cpu_time(t_end)
            write(*,"(12X,A,I4,1X,A,F8.3)") "[Eig] Rank", rank,"Call ZHEEVD finished!  Time cost (s): ", t_end-t_start
        end if

        if (info /= 0) then
            write(*,*) "[Error] Call ZHEEVD failed!"
        end if

        deallocate(work8,rwork8,iwork)

    end subroutine compute_eig_full

    subroutine compute_eig(Hk, eig, jobz)
        use interface, only: zheevr
        use constants, only: rank, prec, timer_debug, band_select, iband
        complex(prec), intent(inout) :: Hk(:,:)
        real(prec), intent(out) :: eig(:)
        character(len=1), intent(in) :: jobz
        
        ! Local variables
        integer :: n, lda, ldz, lwork, lrwork, liwork, info, m
        real(prec) :: vl, vu, abstol
        integer :: il, iu
        complex(8), allocatable :: work8(:), Hk_copy(:,:)
        real(8), allocatable :: rwork8(:), w_full(:)
        integer, allocatable :: iwork(:), isuppz(:)
        complex(8) :: workq(1)
        real(8) :: rworkq(1)
        integer :: iworkq(1)
        real :: t_start, t_end
        character :: range
        integer :: eig_size, m_expected
        
        ! Start timer if debug mode is enabled
        if (timer_debug) call cpu_time(t_start)
        
        ! Get matrix dimension
        n = size(Hk, 1)
        lda = max(1, n)
        ldz = n  ! For zheevr, if eigenvectors are needed, ldz must be at least n
        
        ! Set calculation range based on band_select flag
        if (band_select .and. iband(1) > 0 .and. iband(2) <= n) then
            range = 'I'  ! Select eigenvalues by index
            il = iband(1)
            iu = iband(2)
            m_expected = iu - il + 1
        else
            range = 'A'  ! Compute all eigenvalues
            il = 0
            iu = 0
            m_expected = n
        end if
        
        ! Set tolerance (use default value)
        abstol = 0.0_prec
        
        ! Allocate support array (required by zheevr)
        allocate(isuppz(2 * max(1, m_expected)))
        
        ! Allocate eigenvalue array (need n elements, even for partial calculation)
        allocate(w_full(n))
        
        ! Copy matrix (zheevr destroys input matrix)
        allocate(Hk_copy(n, n))
        Hk_copy = Hk
        
        ! First call: workspace query
        lwork = -1
        lrwork = -1
        liwork = -1
        
        call zheevr(jobz, range, 'U', n, Hk_copy, lda, &
                    vl, vu, il, iu, abstol, &
                    m, w_full, Hk, ldz, isuppz, &
                    workq, lwork, rworkq, lrwork, &
                    iworkq, liwork, info)
        
        ! Check if workspace query succeeded
        if (info /= 0) then
            write(*,*) "[Error] ZHEEVR workspace query failed with info =", info
            deallocate(Hk_copy, isuppz, w_full)
            return
        end if
        
        ! Allocate workspace based on query results
        lwork = max(1, int(real(workq(1))))
        lrwork = max(1, int(rworkq(1)))
        liwork = max(1, iworkq(1))
        
        allocate(work8(lwork), rwork8(lrwork), iwork(liwork))
        
        ! Second call: actual computation
        call zheevr(jobz, range, 'U', n, Hk_copy, lda, &
                    vl, vu, il, iu, abstol, &
                    m, w_full, Hk, ldz, isuppz, &
                    work8, lwork, rwork8, lrwork, &
                    iwork, liwork, info)
        
        ! End timer if debug mode is enabled
        if (timer_debug) then
            call cpu_time(t_end)
            write(*,"(12X,A,I4,1X,A,F8.3)") "[Eig] Rank", rank, &
                  "Call ZHEEVR finished!  Time cost (s): ", t_end - t_start
        end if
        
        ! Check if computation succeeded
        if (timer_debug) then
            if (info /= 0) then
                write(*,*) "[Error] Call ZHEEVR failed with info =", info
            else if (m < m_expected) then
                write(*,"(12X,A,I4,A,I4,A)") "[Warning] ZHEEVR found only ", m, &
                      " eigenvalues, expected ", m_expected, " eigenvalues"
            else
                write(*,"(12X,A,I4,A)") "[Info] ZHEEVR successfully computed ", m, &
                      " eigenvalues"
            end if
        end if
        
        ! Copy results to output array
        if (band_select .and. iband(1) > 0 .and. iband(2) <= n) then
            ! Partial band calculation: copy only computed eigenvalues
            eig_size = min(size(eig), m)
            eig(1:eig_size) = w_full(1:eig_size)
            
            ! Output debug information if needed
            if (timer_debug .and. size(eig) < m) then
                write(*,"(16X,A,I4,A,I4)") "Note: Output array size ", size(eig), &
                      " < computed eigenvalues ", m
            end if
        else
            ! Full band calculation: copy all eigenvalues
            eig_size = min(size(eig), n)
            eig(1:eig_size) = w_full(1:eig_size)
            
            ! Check if output array is large enough
            if (size(eig) < n) then
                write(*,"(12X,A,I4,A,I4)") "[Warning] Output eig array too small: ", &
                      size(eig), " < ", n
            end if
        end if
        
        ! Clean up workspace
        deallocate(work8, rwork8, iwork, isuppz, Hk_copy, w_full)
        
    end subroutine compute_eig

    subroutine calculate_k(kvec,eig,wavef)
        use constants, only: prec,nions
        use tbmodel, only: build_H
        real(prec), intent(in) :: kvec(3)
        real(prec), intent(out) :: eig(:)       
        complex(prec), intent(out), optional :: wavef(:,:)

        complex(prec), allocatable:: Hk(:,:)
        
        if (not(present(wavef))) then
            ! N for eigenvalues only , V also for eigenvectors
            allocate(Hk(nions,nions))
            call build_H(kvec,Hk)
            call compute_eig(Hk, eig, 'N')
            deallocate(Hk)
        else
            call build_H(kvec,wavef)
            call compute_eig(wavef, eig, 'V')
        end if    

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