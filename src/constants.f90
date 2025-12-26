module constants
    implicit none
    logical, public :: timer_mpi
    logical, public :: timer_cpu
    logical, public :: timer_kcalc 
    logical, public :: timer_debug 
    logical, public :: write_hr = .false.
    integer, parameter, public :: prec=SELECTED_REAL_KIND(15)
    real(prec), parameter, public :: pi = 3.14159265358979323846
    real(prec), parameter, public :: k_B = 8.617333262145179e-05 ! Boltzmann constant in eV/K
    real(prec), parameter, public :: epsilon_0_SI = 8.8541878188e-12 ! 
    real(prec), parameter, public :: epsilon_r = 4.9 ! from PhysRevB.102.125403
    real(prec), parameter, public :: elec = 1.602176634e-19
    ! hoppings
    real(prec), parameter, public :: onsite = -0.78_prec ! on-site energy (eV)
    real(prec), parameter, public :: gamma0 = 2.7_prec ! 1-st NN hopping (eV)
    real(prec), parameter, public :: gamma1 = 0.48_prec ! interlayer coupling (eV)
    ! geometric parameters (tbG)
    real(prec), parameter :: a0 = 2.4684713048999796_prec ! lattice constant (AA)
    real(prec), parameter :: a_d = a0/3_prec**0.5_prec ! C-C distance (AA)
    real(prec), parameter :: d0 = 3.39602082978717_prec ! layer distance (AA)
    real(prec), parameter :: l_c = 0.265_prec ! cutoff ratio (AA)
    real(prec), parameter :: r_c = a0*2.5_prec ! layer distance
    ! geometric parameters (plasmon)
    logical :: calc_iqr
    real(prec) :: delta ! broadening
    real(prec) :: efermi = 0.0
    real(prec) :: temperature = 300
    ! decay constants
    real(prec), parameter :: f_decay = log(0.1_prec)/(a_d-a0)
    real(prec), parameter :: q_sigma = d0 * f_decay
    real(prec), parameter :: q_pi = a_d * f_decay

    ! global variants
    logical :: band_select
    logical :: kpt_select
    logical :: job_band
    logical :: write_band
    character(len=3) :: model_type
    character(len=64) :: f_input
    character(len=64) :: f_poscar
    character(len=64) :: f_eig
    character(len=64) :: f_kpoint
    character(len=64) :: f_wannier
    character(len=1) :: q_tag,job_type
    real(prec) :: basis(3,3), basis_rec(3,3)
    real(prec), allocatable :: position_frac(:,:), position_cart(:,:),omegalist(:),qlist(:,:)
    integer :: nk_path 
    integer :: nions,nbands,nomega,iband(2)

    integer ::  rank, nproc, ncache
contains
    subroutine init_constants()
        timer_mpi = .true.
        timer_cpu = .true.
        timer_kcalc = .false.
        timer_debug = .false.
        band_select = .false.
        kpt_select = .false.
        job_band = .false.
        write_band = .false.
        f_input='input.in'
        f_poscar='cont.vasp'
        f_eig='tb_band.dat'
        f_kpoint='KPOINTS'
        nk_path = 10
        ncache = 0
        nomega = 21
        q_tag = "D"
        delta = 1e-3
        calc_iqr = .false.
        model_type = "tbg"
        f_wannier = 'wannier90'
        call set_omegalist(0._prec,0.1_prec)

    end subroutine init_constants

    subroutine set_omegalist(omega_i,omega_f)
        real(prec), intent(in) :: omega_i,omega_f
        integer :: i
        if (allocated(omegalist)) deallocate(omegalist)
        allocate(omegalist(nomega))
        do i=1, nomega
            omegalist(i)=omega_i+(omega_f-omega_i)*(i-1)/(nomega-1)
        end do
    end subroutine set_omegalist

    function get_label(k_frac) result(label)
        real(prec), intent(in) :: k_frac(3)
        real(prec), parameter :: tol = 1.0e-5_prec
        character(len=3) :: label

        if (all(abs(k_frac - [0.0_prec, 0.0_prec, 0.0_prec]) < tol)) then
            label = ' Γ '
        else if (all(abs(k_frac - [1.0_prec/3.0_prec, 1.0_prec/3.0_prec, 0.0_prec]) < tol)) then
            label = ' K '
        else if (all(abs(k_frac - [0.5_prec, 0.0_prec, 0.0_prec]) < tol)) then
            label = ' M '
        else
            label = ' ? '
        end if

    end function get_label 

    subroutine kpath_default(paths)
        real(prec), intent(out), allocatable :: paths(:,:,:)

        if (not(allocated(paths))) then
            allocate(paths(3,2,3))
        end if

        paths(1,1,:) = [0.0_prec, 0.0_prec, 0.0_prec]
        paths(1,2,:) = [1.0_prec/3.0_prec, 1.0_prec/3.0_prec, 0.0_prec]
        paths(2,1,:) = [1.0_prec/3.0_prec, 1.0_prec/3.0_prec, 0.0_prec]
        paths(2,2,:) = [0.5_prec, 0.0_prec, 0.0_prec]
        paths(3,1,:) = [0.5_prec, 0.0_prec, 0.0_prec]
        paths(3,2,:) = [0.0_prec, 0.0_prec, 0.0_prec]

    end subroutine kpath_default
    
end module constants