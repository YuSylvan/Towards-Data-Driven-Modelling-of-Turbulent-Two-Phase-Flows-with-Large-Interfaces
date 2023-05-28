    #if !_MPI
    master = true;
    #else
    master = (mpi_rank == 0);
    #endif
