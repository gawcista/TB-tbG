module tbmodel
    
    implicit none
    
contains
    function f_cutoff(r) result(F_c)
        ! returns the smooth cutoff function
        use constants, only: prec,r_c,l_c
        real(prec), intent(in) :: r
        real(prec) :: F_c
        
        F_c = 1.0_prec/(1.0_prec + exp((r-r_c)/l_c))

    end function f_cutoff

    function tbg_hopping(r_cart) result(t)
        use constants, only: prec, gamma0,gamma1,q_pi,q_sigma,d0,a_d
        real(prec), intent(in) :: r_cart(3)
        real(prec) :: t,V_pi,V_sigma,F_c,n,r

        r = NORM2(r_cart)
        n = r_cart(3)/r
        F_c = f_cutoff(r)
        V_pi = -gamma0*exp(q_pi*(1.0_prec-r/a_d))*F_c
        V_sigma = gamma1*exp(q_sigma*(1.0_prec-r/d0))*F_c
        t = n**2*V_sigma + (1.0_prec-n**2)*V_pi
        
    end function tbg_hopping

    subroutine build_H_tbg(kvec,Hk)
        use constants, only: prec,position_frac,onsite,pi
        use utils, only: r_reduced,r2cart
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
                    HR = tbg_hopping(r2cart(R_vec))   
                    phase = 2.0_prec * pi * dot_product(kvec, R_vec)
                    Hk(i,j) = HR * exp(cmplx(0.0_prec, phase))
                    Hk(j,i) = conjg(Hk(i,j))    
                end if
            end do
        end do  
    end subroutine build_H_tbg

    subroutine build_H(kvec,Hk)
        use constants, only: prec,model_type
        real(prec), intent(in) :: kvec(3)
        complex(prec), intent(out) :: Hk(:,:)

        select case (model_type)
            case('tbg')
                call build_H_tbg(kvec,Hk)
        end select
    end subroutine build_H

    subroutine generate_kpath(paths,nk,k_list_frac,k_list_cart,x_list,labels,lable_positions)  
        use utils, only:k2cart,get_label
        use constants, only: prec
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
    
    
end module tbmodel