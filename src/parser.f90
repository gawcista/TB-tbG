module parser
    use constants
    implicit none
contains
    function isComment(char) result(flag)
        character(len=1) :: char
        logical :: flag

        flag = .false. 
        if (char == '!') flag = .true.
        if (char == '#') flag = .true.
        if (char == '/') flag = .true.

    end function

    subroutine arg_parser()
    end subroutine

    subroutine input_parser(filename)
        use utils, only: generate_kline
        character(len=64), intent(in), optional:: filename
        character(len=64) :: fin,line,key,tmp
        integer :: nqpts,io,unit=101
        real(8) :: omega_i,omega_f,q_start(3),q_end(3)
        
        if (present(filename)) then
            fin = filename
        else
            fin = f_input
        end if

        open(unit=unit, file=fin, status="old")
            do while (.not. eof(unit))
                read(unit,("(A)")) line
                read(line,*,iostat=io) key
                if (.not. isComment(key(1:1))) then
                    select case (key)
                        case ('JOB')
                            read(line,*,iostat=io) tmp,job_type
                        case ('FPOSCAR')
                            read(line,*,iostat=io) tmp,f_poscar
                        case ('FEIGEN')
                            job_band = .true.
                            write_band = .true.
                            read(line,*,iostat=io) tmp,f_eig
                        case ('FKPOINTS')
                            kpt_select = .true.
                            read(line,*,iostat=io) tmp,f_kpoint
                        case ('FWANNIER')
                            model_type = 'wan'
                            read(line,*,iostat=io) tmp,f_wannier
                        case ('IBAND')
                            band_select = .true.
                            read(line,*,iostat=io) tmp,iband
                            nbands = iband(2)-iband(1)+1
                        case ('NCACHE')
                            read(line,*,iostat=io) tmp,ncache
                        case ('OMEGA')
                            read(line,*,iostat=io) tmp,nomega,omega_i,omega_f
                            call set_omegalist(omega_i,omega_f)
                        case ('QPOINT')
                            read(line,*,iostat=io) tmp,nqpts,q_start,q_end
                            call generate_kline(nqpts,q_start,q_end,qlist)
                        case ('QTAG')
                            read(line,*,iostat=io) tmp,q_tag
                        case ('EFERMI')
                            read(line,*,iostat=io) tmp,efermi
                        case ('TEMPERATURE')
                            read(line,*,iostat=io) tmp,temperature
                        case ('DELTA')
                            read(line,*,iostat=io) tmp,delta
                        case('CALC_IQR')
                            read(line,*,iostat=io) tmp,calc_iqr
                        case('TMPI')
                            read(line,*,iostat=io) tmp,timer_mpi
                        case('TDEBUG')
                            read(line,*,iostat=io) tmp,timer_debug
                        case('MODEL')
                            read(line,*,iostat=io) tmp,model_type         
                    end select                        
                                 
                end if
            end do

    

        close(unit)

        

    end subroutine
end module parser