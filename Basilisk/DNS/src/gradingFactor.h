double gradingFactor
(
    const double startWidth,
    const double length,
    const int M
)
{
    // Function to compute the geometric grading factor given a startwidth of
    // the cell, an interval length and M number of cells.

    const double Niter = 1024;
    const double TOL = 1e-10;

    // Solve in log-space using the transformation y = log(f-1). This is
    // numerically better conditioned. The initial guess is f = 1.1.

    double f = 1.1;
    double y = log(f-1.0);

    for (int i = 0; i < Niter; i++)
    {
        // Update increment dy = g(y)/(dg(y)/dy)

        const double dy =
            (exp(-y)*(pow(f,M) - 1.0) - length/startWidth)
          / (M*pow(f,M-1) - exp(-y)*(pow(f,M)-1.0));

        y -= dy;

        f = exp(y)+1.0;

        if (fabs(dy) < TOL)
        {
            break;
        }

        if (i == Niter-1)
        {
            #if _MPI
            if (mpi_rank == 0)
            #endif
            {
                fprintf(stdout, "Could not find grading factor f\n");
            }

            exit(1);
        }
    }

    return f;
}
