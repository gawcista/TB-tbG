module eels 
    implicit none
contains
    function Fermi_Dirac(E,Ef,T) result(f)
        use constants, only: prec, efermi, temperature ,k_B
        real(prec):: E,f
        real(prec), optional :: Ef,T
        real (prec) :: E_fermi
        real (prec) :: temp
        
        E_fermi = efermi
        temp = temperature
        if (present(Ef)) E_fermi = Ef
        if (present(T)) temp = T
        f = 1._prec/(exp((E-E_fermi)/(k_B*temp)) + 1._prec)

    end function Fermi_Dirac

    function V_q(q) result(V)
        ! return e²/((2pi)^2*ϵ₀*ϵ_r*q)
        use constants, only:elec,epsilon_0_SI,prec,epsilon_r,pi
        real(prec) :: q(3),V,V_f,q_mod
        V_f = elec*1.e10/epsilon_0_SI ! e²/ϵ₀~45.2 eV*angstrom
        q_mod = norm2(q)
        if (q_mod<1e-5) q_mod = 1e-5
        V = V_f/((2*pi)**2*epsilon_r*q_mod) 
    end function V_q 

    subroutine calculate_IPF_k(omegalist,eig_k,eig_kq,phi_k,phi_kq,resultlist,q_cart,matrix_contrib)
        ! require cartesian q
        use constants, only: prec,nbands,nions,write_matrix
        real(prec), intent(in) :: omegalist(:),eig_k(:),eig_kq(:)
        complex(prec), intent(in) :: phi_k(:,:),phi_kq(:,:)
        real(prec), intent(in), optional :: q_cart(3)
        complex(prec), intent(inout) :: resultlist(:)
        complex(prec), intent(inout), optional :: matrix_contrib(:,:,:)  ! (nbands,nbands,nomega)
        complex(prec), allocatable :: phase(:),phi_k_scaled(:,:),M(:,:)
        real(prec), allocatable :: f_k(:), f_kq(:)
        integer :: i

        allocate(phase(nions),phi_k_scaled(nions,nbands),M(nbands,nbands))
        allocate(f_k(nbands),f_kq(nbands))
        do i = 1, nbands
            f_k(i) = Fermi_Dirac(eig_k(i))
            f_kq(i) = Fermi_Dirac(eig_kq(i))
        end do

        if (present(matrix_contrib) .and. write_matrix) matrix_contrib = (0.0_prec, 0.0_prec)

        if (present(q_cart)) then
            call calculate_IPF_k_work(omegalist,eig_k,eig_kq,f_k,f_kq,phi_k,phi_kq,resultlist, &
                                      phase,phi_k_scaled,M,q_cart,matrix_contrib)
        else if (present(matrix_contrib)) then
            call calculate_IPF_k_work(omegalist,eig_k,eig_kq,f_k,f_kq,phi_k,phi_kq,resultlist, &
                                      phase,phi_k_scaled,M,matrix_contrib=matrix_contrib)
        else
            call calculate_IPF_k_work(omegalist,eig_k,eig_kq,f_k,f_kq,phi_k,phi_kq,resultlist, &
                                      phase,phi_k_scaled,M)
        end if

        deallocate(phase,phi_k_scaled,M,f_k,f_kq)

    end subroutine calculate_IPF_k

    subroutine calculate_IPF_k_work(omegalist,eig_k,eig_kq,f_k,f_kq,phi_k,phi_kq,resultlist, &
                                   phase,phi_k_scaled,M,q_cart,matrix_contrib)
        ! require cartesian q
        use constants, only: prec,nbands,nions,position_cart,delta,eels_mode,write_matrix
        use interface, only: zgemm
        real(prec), intent(in) :: omegalist(:),eig_k(:),eig_kq(:),f_k(:),f_kq(:)
        complex(prec), intent(in) :: phi_k(:,:),phi_kq(:,:)
        complex(prec), intent(inout) :: resultlist(:)
        complex(prec), intent(inout) :: phase(:),phi_k_scaled(:,:),M(:,:)
        real(prec), intent(in), optional :: q_cart(3)
        complex(prec), intent(inout), optional :: matrix_contrib(:,:,:)  ! (nbands,nbands,nomega)
        integer :: i,iband_k,iband_kq,nomega
        real(prec) :: delta_f, delta_E, matrix_weight

        if (present(q_cart)) then
            do i = 1,nions
                phase(i) = exp(cmplx(0.0_prec, dot_product(q_cart,position_cart(i,:)),prec))
            end do
        else
            phase = cmplx(1.0_prec, 0.0_prec, prec)
        end if

        do i=1, nbands
            phi_k_scaled(:,i) = phi_k(:,i) * phase
        end do

        call zgemm('C', 'N', nbands, nbands, nions, &
               cmplx(1.0_prec, 0.0_prec, prec), &
               phi_kq, nions, &
               phi_k_scaled, nions, &
               cmplx(0.0_prec, 0.0_prec, prec), &
               M, nbands)

        nomega = size(omegalist)

        !$OMP PARALLEL DO PRIVATE(delta_f, delta_E, matrix_weight, i) REDUCTION(+:resultlist) &
        !$OMP& COLLAPSE(2) SCHEDULE(static)
        do iband_k=1,nbands
            do iband_kq=1,nbands
                if (eels_mode==1 .and. iband_k/=iband_kq) cycle
                if (eels_mode==2 .and. iband_k==iband_kq) cycle
                delta_f = f_k(iband_k) - f_kq(iband_kq)
                delta_E = eig_k(iband_k) - eig_kq(iband_kq)
                matrix_weight = real(conjg(M(iband_kq, iband_k)) * M(iband_kq, iband_k), prec)
                do i = 1, nomega
                resultlist(i) = resultlist(i) + matrix_weight * delta_f / &
                                (cmplx(delta_E + omegalist(i),delta,prec))
                end do
                ! Store individual transition contributions if needed
                if (present(matrix_contrib) .and. write_matrix) then
                    do i = 1, nomega
                        matrix_contrib(iband_kq, iband_k, i) = &
                            matrix_contrib(iband_kq, iband_k, i) + &
                            matrix_weight * delta_f / (cmplx(delta_E + omegalist(i),delta,prec))
                    end do
                end if
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine calculate_IPF_k_work

    subroutine calculate_IPF_klist(omegalist,klist_frac,q_list_frac,IPF,matrix_contrib)
        use constants, only: prec,calc_iqr,write_matrix,nbands,nions
        use utils, only: k2frac,k2cart,k2frac_list
        use solver, only : calculate_klist
        real(prec), intent(in) :: omegalist(:),klist_frac(:,:),q_list_frac(:,:)
        complex(prec), allocatable, intent(inout):: IPF(:,:)
        complex(prec), allocatable, intent(inout), optional :: matrix_contrib(:,:,:,:)  ! (nbands,nbands,nomega,nqpts)
        real(prec) :: q_frac(3),q_cart(3)
        real(prec), allocatable :: f_k(:), f_kq(:)
        real(prec), allocatable :: kqlist_frac(:,:),eig_k(:,:),eig_kq(:,:)
        complex(prec),allocatable :: phi_k(:,:,:),phi_kq(:,:,:),phase(:),phi_k_scaled(:,:),M(:,:)
        integer :: nkpts,nqpts,ik,iq,nomega,iband

        nkpts = size(klist_frac,1)
        nqpts = size(q_list_frac,1)
        nomega = size(omegalist)
        allocate(kqlist_frac(nkpts,3))
        allocate(phase(nions),phi_k_scaled(nions,nbands),M(nbands,nbands))
        allocate(f_k(nbands),f_kq(nbands))
        if (.not. allocated(IPF)) then
            allocate(IPF(nomega,nqpts))
        end if 
        IPF = (0.0_prec, 0.0_prec)
        ! Allocate matrix_contrib if needed
        if (present(matrix_contrib) .and. write_matrix) then
            if (.not. allocated(matrix_contrib)) then
                allocate(matrix_contrib(nbands, nbands, nomega, nqpts))
            end if
            matrix_contrib = (0.0_prec, 0.0_prec)
        end if
        call calculate_klist(klist_frac,eig_k,phi_k)
        do iq = 1, nqpts
            q_frac = q_list_frac(iq,:)
            q_cart = k2cart(q_frac)
            do ik = 1, nkpts
                kqlist_frac(ik,:) = klist_frac(ik,:) + q_frac
            end do
            call calculate_klist(kqlist_frac,eig_kq,phi_kq)
            do ik=1,nkpts
                do iband = 1, nbands
                    f_k(iband) = Fermi_Dirac(eig_k(iband,ik))
                    f_kq(iband) = Fermi_Dirac(eig_kq(iband,ik))
                end do
                if (present(matrix_contrib) .and. write_matrix) then
                    if (calc_iqr) then
                        call calculate_IPF_k_work(omegalist,eig_k(:,ik),eig_kq(:,ik),f_k,f_kq, &
                                                  phi_k(:,:,ik),phi_kq(:,:,ik),IPF(:,iq), &
                                                  phase,phi_k_scaled,M,q_cart,matrix_contrib(:,:,:,iq))
                    else
                        call calculate_IPF_k_work(omegalist,eig_k(:,ik),eig_kq(:,ik),f_k,f_kq, &
                                                  phi_k(:,:,ik),phi_kq(:,:,ik),IPF(:,iq), &
                                                  phase,phi_k_scaled,M,matrix_contrib=matrix_contrib(:,:,:,iq))
                    end if
                else
                    if (calc_iqr) then
                        call calculate_IPF_k_work(omegalist,eig_k(:,ik),eig_kq(:,ik),f_k,f_kq, &
                                                  phi_k(:,:,ik),phi_kq(:,:,ik),IPF(:,iq), &
                                                  phase,phi_k_scaled,M,q_cart)
                    else
                        call calculate_IPF_k_work(omegalist,eig_k(:,ik),eig_kq(:,ik),f_k,f_kq, &
                                                  phi_k(:,:,ik),phi_kq(:,:,ik),IPF(:,iq), &
                                                  phase,phi_k_scaled,M)
                    end if
                end if
            end do
            deallocate(eig_kq,phi_kq)
        end do
        deallocate(eig_k,phi_k,kqlist_frac,phase,phi_k_scaled,M,f_k,f_kq)

    end subroutine calculate_IPF_klist

    subroutine calculate_eels(IPF,q_list_cart,resultlist)
        use constants, only : prec
        complex(prec), intent(in) :: IPF(:,:)
        real(prec), intent(in) :: q_list_cart(:,:)
        real(prec), allocatable, intent(out) :: resultlist(:,:)
        complex(prec) :: epsilon
        integer :: nomega,nqpts,iomega,iq

        nomega = size(IPF,1)
        nqpts = size(IPF,2)
        if (.not. allocated(resultlist)) allocate(resultlist(nomega,nqpts))
        do iq = 1, nqpts
            do iomega = 1, nomega
                epsilon = 1._prec - V_q(q_list_cart(iq,:)) * IPF(iomega,iq)
                resultlist(iomega,iq) = - imag(1._prec/epsilon)
            end do
        end do

    end subroutine calculate_eels

end module eels
