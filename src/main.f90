program tbgsolver
    use constants, only: job_type,qlist,q_tag,rank,nbands
    use mpi_f08
    use mpi_solver
    implicit none
    integer :: ierr


    call init_mpi()
    call init_POSCAR_mpi()
    select case(job_type)
    case("B")
        if (rank==0) write(*,"(A,I4,1X,A)") "[Main] Job type set to band structure, ",nbands,"bands would be calculated."
        call calculate_band_mpi()
    case("E")
        if (rank==0) write(*,"(A)") "[Main] Job type set to EELS"
        call calculate_eels_mpi(qlist,q_tag)
    end select
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)

end program tbgsolver