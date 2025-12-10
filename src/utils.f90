module utils
    use constants
    implicit none

contains
    function cross2d(a,b) result(c)
        real(prec), intent(in) :: a(2),b(2)
        real(prec) :: c
        c = a(1) * b(2) - a(2) * b(1)     
    end function cross2d

    function cross3d(a, b) result(c)
        implicit none
        real(prec), intent(in) :: a(3), b(3)
        real(prec) :: c(3)
  
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3d

    subroutine matrix_inv(A, A_inv)
        use interface, only:dgetrf,dgetri
        real(prec), intent(in) :: A(:,:)
        real(prec), intent(out) :: A_inv(:,:)
        integer :: n, ipiv(size(A,1)), info
        real(prec) :: work(size(A,1))
        
        n = size(A,1)
        A_inv = A
        call dgetrf(n, n, A_inv, n, ipiv, info)
        call dgetri(n, A_inv, n, ipiv, work, n, info)
        
    end subroutine matrix_inv
    
    function count_words(line) result(n)
        character(len=*), intent(in) :: line
        integer :: n
        integer :: i, length
        logical :: in_word

        n = 0
        in_word = .false.
        length = len_trim(line)

        do i = 1, length
            if (line(i:i) /= ' ') then
                if (.not. in_word) then
                    n = n + 1
                    in_word = .true.
                end if
            else
                in_word = .false.
            end if
        end do
    end function count_words

    subroutine get_reciprocal_basis(basis_r, basis_k)
        real(prec), intent(in) :: basis_r(3,3)
        real(prec), intent(out) :: basis_k(3,3)
        real(prec) :: volume

        ! cell volume = a1 · (a2 × a3)
        volume = dot_product(basis_r(1,:), cross3d(basis_r(2,:), basis_r(3,:)))
        ! reciprocal basis vectors
        basis_k(1,:) = 2*pi * cross3d(basis_r(2,:), basis_r(3,:)) / volume
        basis_k(2,:) = 2*pi * cross3d(basis_r(3,:), basis_r(1,:)) / volume  
        basis_k(3,:) = 2*pi * cross3d(basis_r(1,:), basis_r(2,:)) / volume
        
    end subroutine get_reciprocal_basis

    function r_reduced(r_frac) result(r_R)
        real(prec), intent(in) :: r_frac(3)
        real(prec) :: r_R(3)
        integer :: i

        r_R = r_frac
        do i = 1,2
            if (r_R(i) >= 0.5_prec) then
                r_R(i) = r_R(i) - 1.0_prec
            else if (r_R(i) < -0.5_prec) then
                r_R(i) = r_R(i) + 1.0_prec
            end if
        end do

    end function r_reduced

    function r2cart_base(basis_r,pos_frac) result(pos_cart)
        real(prec), intent(in) :: basis_r(3,3), pos_frac(3)
        real(prec) :: pos_cart(3)

        pos_cart = matmul(transpose(basis_r), pos_frac)

    end function r2cart_base

    function r2frac_base(basis_r,pos_cart) result(pos_frac)
        real(prec), intent(in) :: basis_r(3,3), pos_cart(3)
        real(prec) :: pos_frac(3)
        real(prec) :: inv_basis(3,3)

        call matrix_inv(transpose(basis_r), inv_basis)
        pos_frac = matmul(inv_basis, pos_cart)
        
    end function r2frac_base

    function r2cart(pos_frac) result(pos_cart)
        real(prec), intent(in) :: pos_frac(3)
        real(prec) :: pos_cart(3)

        pos_cart = r2cart_base(basis, pos_frac)

    end function r2cart

    function r2frac(pos_cart) result(pos_frac)
        real(prec), intent(in) :: pos_cart(3)
        real(prec) :: pos_frac(3)

        pos_frac = r2frac_base(basis,pos_cart)
        
    end function r2frac
    

    function k2cart(k_frac) result(k_cart)
        real(prec), intent(in) :: k_frac(3)
        real(prec) :: k_cart(3)

        k_cart = matmul(k_frac, basis_rec)
        
    end function k2cart

    function k2frac(k_cart) result(k_frac)
        real(prec), intent(in) :: k_cart(3)
        real(prec) :: k_frac(3)
        real(prec) :: inv_recip_basis(3,3)

        call matrix_inv(basis_rec, inv_recip_basis)
        k_frac = matmul(k_cart, transpose(inv_recip_basis))
        
    end function k2frac

    subroutine k2cart_list(k_list_frac,k_list_cart)
        real(prec), intent(in) :: k_list_frac(:,:)
        real(prec), intent(out) :: k_list_cart(size(k_list_frac,1),3)
        integer :: i, nkpoints

        nkpoints = size(k_list_frac,1)
        do i = 1,nkpoints
            k_list_cart(i,:) = k2cart(k_list_frac(i,:))
        end do

    end subroutine k2cart_list

    subroutine k2frac_list(k_list_cart,k_list_frac)
        real(prec), intent(in) :: k_list_cart(:,:)
        real(prec), intent(out) :: k_list_frac(size(k_list_cart,1),3)
        integer :: i, nkpoints

        nkpoints = size(k_list_cart,1)
        do i = 1,nkpoints
            k_list_frac(i,:) = k2frac(k_list_cart(i,:))
        end do

    end subroutine k2frac_list

    subroutine generate_mesh_gamma(nk,klist_frac)
        integer, intent(in) :: nk(3)
        real(prec), allocatable, intent(out) :: klist_frac(:,:)
        real(prec) :: kx,ky,kz
        integer :: i, j, k, idx
        

        if (.not. allocated(klist_frac)) allocate(klist_frac(nk(1)*nk(2)*nk(3), 3))
        idx = 0

        do i = 1, nk(1)
            kx = (i - 1) / real(nk(1), prec)
            do j = 1, nk(2)
                ky = (j - 1) / real(nk(2), prec)
                do k = 1, nk(3)
                    kz = (k - 1) / real(nk(3), prec)
                    idx = idx + 1
                    klist_frac(idx, :) = [kx, ky, kz]
                end do
            end do
        end do
    end subroutine generate_mesh_gamma

    subroutine generate_mesh_mp(nk,klist_frac)
        integer, intent(in) :: nk(3)
        real(prec), allocatable, intent(out) :: klist_frac(:,:)
        real(prec) :: kx,ky,kz
        integer :: i, j, k, idx
        

        if (.not. allocated(klist_frac)) allocate(klist_frac(nk(1)*nk(2)*nk(3), 3))
        idx = 0
        
        do i = 1, nk(1)
            kx = (2*i - nk(1) - 1) / (2.0_prec * nk(1))
            do j = 1, nk(2)
                ky = (2*j - nk(2) - 1) / (2.0_prec * nk(2))
                do k = 1, nk(3)
                    kz = (2*k - nk(3) - 1) / (2.0_prec * nk(3))
                    idx = idx + 1
                    klist_frac(idx, :) = [kx, ky, kz]
                end do
            end do
        end do
    end subroutine generate_mesh_mp

    subroutine generate_kline(nk,k_start,k_end,klist)
        integer, intent(in) :: nk
        real(prec), intent(in) :: k_start(3),k_end(3)
        real(prec), intent(out), allocatable :: klist(:,:)
        integer :: i
        if (allocated(klist)) deallocate(klist)
        allocate(klist(nk,3))
        if (nk > 1) then
            do i=1,nk
                klist(i,:) = k_start + (i-1._prec)*(k_end-k_start)/(nk-1)
            end do 
        else
            klist(1,:) = k_start
        end if
    end subroutine generate_kline

    

    
end module utils