module ioutils
    use utils
    implicit none
    
contains
    subroutine readPOSCAR(filename,lattice)
        use constants, only:prec
        use types, only: crystal
        character(len=64), intent(in):: filename
        type(crystal), intent(out):: lattice

        integer :: i,nelements
        character(len=64) :: temp_line
        character(len=1) :: tag
        real(prec) :: temp_position(3)
        real :: t_start,t_end
        integer :: unit=100

        if (timer_cpu) then
            write(*, '(A,1X,A,1X,A,1X,$)') "[IO] Reading POSCAR file: ", trim(filename), "..."
            call cpu_time(t_start)
        end if

        open(unit, file=filename, status='old', action='read')
            read(unit, '(A)') lattice%title
            read(unit, *) lattice%frac
            do i = 1, 3
                read(unit, *) lattice%basis(i,:)
            end do
            lattice%basis = lattice%frac * lattice%basis
            call get_reciprocal_basis(lattice%basis, lattice%basis_rec)
            read(unit, '(A)') temp_line ! read elements line
            nelements = count_words(temp_line)
            allocate(lattice%elements(nelements))
            allocate(lattice%natoms(nelements))
            read(temp_line, *) lattice%elements
            read(unit, *) lattice%natoms
            lattice%nions = sum(lattice%natoms)

            allocate(lattice%atom(lattice%nions))
            read(unit, '(A)') temp_line ! read selective dynamics or coordinate type
            tag = adjustl(temp_line)
            if (tag(1:1)=='S') then
                lattice%selective_dynamics = .true.
                read(unit, '(A)') temp_line
                tag = trim(adjustl(temp_line))
            else 
                lattice%selective_dynamics = .false.
            end if
            lattice%tag = tag(1:1)           
            do i = 1, lattice%nions
                if (lattice%selective_dynamics) then
                    read(unit,*) temp_position,lattice%atom(i)%dynamics
                else
                    read(unit,*) temp_position
                end if
                select case(lattice%tag)
                    case('D')
                        lattice%atom(i)%position_frac = temp_position
                        lattice%atom(i)%position_cart = r2cart_base(lattice%basis,temp_position)
                    case('C')
                        lattice%atom(i)%position_cart = temp_position
                        lattice%atom(i)%position_frac = r2frac_base(lattice%basis,temp_position)
                end select
            end do
        close(unit)
        
        if (timer_cpu) then
            call cpu_time(t_end)
            write(*, '(A,1X,F8.3)') "Done! Time elapsed (s): ", t_end - t_start
        end if

    end subroutine readPOSCAR

    subroutine parseCrystal(lattice)
        use types, only: crystal,nbands,band_select
        use constants, only: position_cart,position_frac,nions,basis,basis_rec
        type(crystal), intent(in) :: lattice
        integer :: i

        nions = lattice%nions
        basis = lattice%basis
        if (.not. band_select) then
            iband(1) = 1
            iband(2) = nions
        end if
        nbands = iband(2) - iband(1) + 1
        basis_rec = lattice%basis_rec
        allocate(position_frac(nions,3))
        allocate(position_cart(nions,3))
        do i = 1,nions
            position_frac(i,:) = lattice%atom(i)%position_frac
            position_cart(i,:) = lattice%atom(i)%position_cart
        end do

    end subroutine parseCrystal

    subroutine parsePOSCAR(filename)
        use types, only: crystal
        character(len=64), intent(in) :: filename
        type(crystal) :: lattice

        call readPOSCAR(filename,lattice)
        call parseCrystal(lattice)

    end subroutine parsePOSCAR

    subroutine readKPOINTS(filename,klist_frac)
        use constants, only:prec,f_kpoint
        type :: node
            real(prec) :: data(3)
            type(node), pointer :: next => null()
        end type node
        character(len=64), optional, intent(in):: filename
        real(prec), allocatable,optional, intent(out) :: klist_frac(:,:)
        character(len=64) :: fin,mode,tag
        character(len=128) :: line
        integer :: i,unit,nk(3),nkpts
        type(node), pointer :: head, current
        real(prec) :: values(3)

        unit = 102

        if (present(filename)) then
            fin = filename
        else 
            fin = f_kpoint
        end if

        nkpts = 0
        open(unit, file=fin, status='old', action='read')
            read(unit, '(A)') line
            read(unit, '(I)') nk(1)
            if (nk(1) == 0) then
                ! determine k-mesh automatically
                read(unit, '(A)') mode
                read(unit, *) nk
                if (mode(1:1)=='G' .or. mode(1:1)=='g') then
                    call generate_mesh_gamma(nk,klist_frac)
                else if (mode(1:1)=='M' .or. mode(1:1)=='m') then
                    call generate_mesh_mp(nk,klist_frac)
                end if
            else 
                read(unit, '(A)') mode
                if (mode(1:1)=='L' .or. mode(1:1)=='l') then
                    nk_path = nk(1)
                    read(unit, '(A)') tag
                    nullify(head)
                    do while (.not. eof(unit))
                        read(unit, '(A)') line
                        if (len_trim(line)==0) cycle ! skip empty lines
                        read(line,*) values
                        allocate(current)
                        current%data = values
                        current%next => head
                        head => current
                        nkpts = nkpts + 1
                        
                    end do
                end if
            end if 
        close(unit)

        if (nkpts>0) then
            allocate(klist_frac(nkpts, 3))
            current => head
            do i = nkpts, 1, -1
                klist_frac(i, :) = current%data
                current => current%next
                
            end do
            do while (associated(head))
                current => head
                head => head%next
                deallocate(current)
            end do
        end if

    end subroutine readKPOINTS

    subroutine readWannier(seedname)
        character(len=64), optional, intent(in):: seedname
        character(len=64) :: fin
        integer :: unit
        
        unit = 103
        if (present(seedname)) then
            fin = seedname
        else 
            fin = f_wannier
        end if



    end subroutine readWannier


    subroutine writeBand(x_list,eig,style,istart,iend)
        use constants, only:prec,f_eig,nbands
        real(prec), intent(in) :: x_list(:),eig(:,:)
        character(len=3), intent(in), optional :: style
        integer, intent(in), optional :: istart,iend
        character(len=3) :: dstyle = 'gnu'
        character(len=64) :: fout
        integer :: unit, i, j, nkpoints, ib_start, ib_end
        real :: t_start,t_end
        
        call cpu_time(t_start)

        unit = 200
        nkpoints = size(x_list,1)

        if (present(style)) then
            dstyle = style
        end if
        if (present(istart)) then
            ib_start = istart
        else
            ib_start = 1
        end if
        if (present(iend)) then
            ib_end = iend
        else
            ib_end = nbands
        end if
        write(*,'(A,I4,A,I4,1X,A,1X,I5,1X,A,1X,A,1X)') "[IO] Writing band ",ib_start,"-",ib_end,"on",nkpoints,"k-points to",fout
        open(unit, file=f_eig, status='replace', action='write')
            if (dstyle == 'gnu') then
                write(*,'(A,1X,A,1X,$)') "[IO] Output style: gnuplot"
                do i = ib_start,ib_end
                    do j = 1, nkpoints
                        write(unit,'(F12.8,1X,F20.15)') x_list(j), eig(i,j)
                    end do
                    write(unit,*) ! blank line between bands
                end do
            end if
        close(unit)

        call cpu_time(t_end)
        write(*, '(A,1X,F8.3)') "[IO] Done! Time elapsed (s): ", t_end - t_start
    end subroutine writeBand

    subroutine writeLabels(labels, lable_positions,filename)
        character(len=3), intent(in), allocatable :: labels(:)
        real(prec), intent(in), allocatable:: lable_positions(:)
        character(len=64), intent(in), optional :: filename
        character(len=64) :: fout = 'tb_labels.dat'
        integer :: unit,i,n

        unit=300
        n = size(labels,1)
        if (present(filename)) fout = filename
        open(unit, file=fout, status='replace', action='write')
            do i = 1,n
                write(unit,"(A,1X,F12.8)") labels(i),lable_positions(i)
            end do
            
        close(unit)
    end subroutine writeLabels

    subroutine writeEELS(q_list_cart,omega,eloss,filename)
        use constants, only: nomega,temperature,delta,efermi
        real(prec), intent(in) :: q_list_cart(:,:),omega(:),eloss(:,:)
        character(len=64), intent(in), optional :: filename
        character(len=64) :: fout = 'tb_eels.dat'
        integer :: unit,iomega,iq,nqpts

        unit = 400
        nqpts = size(q_list_cart,1)
        if (present(filename)) fout = filename
        write(*,"(A,1X,A)") "[IO] Wrinting EELS data to file:",fout
        
        open(unit, file=fout, status='replace', action='write')
            write(unit,"(A,F8.4,A,F8.4,A,F8.4,A,I8)") "# T=",temperature,"; delta=",delta*1000," meV; E-fermi=",efermi," eV"
            do iq = 1,nqpts
                do iomega = 1,nomega
                    write(unit,'(5F12.8)') q_list_cart(iq,:),omega(iomega),eloss(iomega,iq)
                end do
                write(unit,*)
            end do
        close(unit)

    end subroutine writeEELS

end module ioutils