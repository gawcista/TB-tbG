module tbmodel
    use constants
    use interface, only: zheevd
    use utils

    implicit none
    
contains
    function f_cutoff(r) result(F_c)
        ! returns the smooth cutoff function
        real(prec), intent(in) :: r
        real(prec) :: F_c
        
        F_c = 1.0_prec/(1.0_prec + exp((r-r_c)/l_c))

    end function f_cutoff

    function hopping(r_cart) result(t)
        real(prec), intent(in) :: r_cart(3)
        real(prec) :: t,V_pi,V_sigma,F_c,n,r

        r = NORM2(r_cart)
        n = r_cart(3)/r
        F_c = f_cutoff(r)
        V_pi = -gamma0*exp(q_pi*(1.0_prec-r/a_d))*F_c
        V_sigma = gamma1*exp(q_sigma*(1.0_prec-r/d0))*F_c
        t = n**2*V_sigma + (1.0_prec-n**2)*V_pi
        
    end function hopping

    subroutine build_H(kvec,Hk)
        real(prec), intent(in) :: kvec(3)
        complex(prec), intent(out) :: Hk(:,:)
        integer :: i,j
        real(prec) :: R_vec(3),phase
        complex(prec) :: HR


        do i = 1,size(position_frac,1)
            do j = i,size(position_frac,1)
                if (i==j) then
                    Hk(i,i) = onsite 
                else
                    R_vec = r_reduced(position_frac(j,:) - position_frac(i,:))
                    HR = hopping(r2cart(R_vec))   
                    phase = 2.0_prec * pi * dot_product(kvec, R_vec)
                    Hk(i,j) = HR * exp(cmplx(0.0_prec, phase))
                    Hk(j,i) = conjg(Hk(i,j))    
                end if
            end do
        end do    

    end subroutine build_H

    subroutine compute_eig(Hk, eig, jobz)
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
            write(*,"(12X,A,1X,F8.3)") "[eig] Call ZHEEVD finished!  Time cost (s): ", t_end-t_start
        end if

        if (info /= 0) then
            write(*,*) "[Error] Call ZHEEVD failed!"
        end if

        deallocate(work8,rwork8,iwork)

    end subroutine compute_eig

    subroutine generate_kpath(paths,nk,k_list_frac,k_list_cart,x_list,labels,lable_positions)  
        real(prec), intent(in) :: paths(:,:,:)
        integer, intent(in) :: nk
        real(prec), allocatable, intent(out) :: k_list_frac(:,:),k_list_cart(:,:),x_list(:)
        character(len=3), allocatable, intent(out) :: labels(:)
        real(prec), allocatable, intent(out) :: lable_positions(:)
        integer :: npaths,i,j,nkpts
        real(prec) :: cumulative_distance,k_frac(3),k_cart(3)

        write(*,"(A,1X,$)") "[Main] Generating k-point path ..."

        cumulative_distance = 0.0_prec
        npaths = size(paths,1)
        if ( .not. allocated(k_list_frac) ) then
            allocate(k_list_frac(nk*npaths+1,3))
        end if
        if ( .not. allocated(k_list_cart) ) then
            allocate(k_list_cart(nk*npaths+1,3))    
        end if
        allocate(x_list(nk*npaths+1))
        allocate(labels(npaths+1))
        allocate(lable_positions(npaths+1))
        k_list_frac(1,:) = paths(1,1,:)
        k_list_cart(1,:) = k2cart(k_list_frac(1,:))
        x_list(1) = 0.0_prec
        labels(1) = get_label(paths(1,1,:))
        lable_positions(1) = 0.0_prec
        do i = 1,npaths
            do j = 1,nk
                k_frac = paths(i,1,:) + &
                    (real(j,prec)/real(nk,prec)) * (paths(i,2,:) - paths(i,1,:))
                k_cart = k2cart( k_frac)
                k_list_frac((i-1)*nk + j + 1, :) = k_frac
                k_list_cart((i-1)*nk + j + 1, :) = k_cart
                cumulative_distance = cumulative_distance + &
                    norm2(k_cart - k_list_cart((i-1)*nk + j, :))
                x_list((i-1)*nk + j + 1) = cumulative_distance
            end do
            labels(i+1) = get_label(paths(i,2,:))
            lable_positions(i+1) = cumulative_distance
        end do 
        nkpts = size(k_list_frac,1)
        write(*,"(I4,1X,A)") nkpts,"k-points generated."
    
    end subroutine generate_kpath

    subroutine calculate_k(kvec,eig,wavef)
        real(prec), intent(in) :: kvec(3)
        real(prec), intent(out) :: eig(:)       
        complex(prec), intent(out), optional :: wavef(:,:)

        complex(prec), allocatable:: Hk(:,:)
        
        if (not(present(wavef))) then
            ! N for eigenvalues only , V also for eigenvectors
            allocate(Hk(size(eig,1),size(eig,1)))
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

        if (.not. allocated(eig)) allocate(eig(nions,nkpts))
        if (present(wavef) .and. .not. allocated(wavef)) allocate(wavef(nions,nions,nkpts))

        do i = 1,nkpts
            if (present(wavef)) then
                call calculate_k(kpoints(i,:),eig(:,i),wavef(:,:,i))
            else
                call calculate_k(kpoints(i,:),eig(:,i))
            end if
        end do

    end subroutine calculate_klist
    
    
end module tbmodel