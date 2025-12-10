module types
    use constants
    implicit none
    type atom_unit
        character(len=2) :: element
        real(prec) :: position_frac(3),position_cart(3)
        logical :: dynamics(3)=[.true.,.true.,.true.]
        !character(len=1) :: dynamics(3)=['T','T','T']
    end type atom_unit
    
    type crystal
        character(len=64) :: title
        real(prec) :: basis(3,3),basis_rec(3,3)
        type(atom_unit), allocatable :: atom(:)
        real :: frac = 1.0
        integer, allocatable :: natoms(:)
        character(len=2), allocatable :: elements(:)
        logical :: selective_dynamics 
        integer :: nions=0
        character(len=1) tag
    end type crystal
contains
end module types