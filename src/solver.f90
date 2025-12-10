module solver
    use constants
    use ioutils
    use tbmodel
    implicit none

contains

    subroutine calc_band_default()

        real(prec), allocatable :: eig(:,:)
        real(prec), allocatable:: paths(:,:,:),klist_frac(:,:),klist_cart(:,:),xlist(:),lable_positions(:)
        character(len=3), allocatable :: labels(:)

        call kpath_default(paths)
        call generate_kpath(paths, 10, klist_frac, klist_cart, xlist, labels, lable_positions)
        call calculate_klist(klist_frac,eig)
        call writeBand(xlist,eig)
        call writeLabels(labels, lable_positions)
        
    end subroutine calc_band_default

    !subroutine calc_eels_default()
    !end subroutine calc_eels_default()
end module solver