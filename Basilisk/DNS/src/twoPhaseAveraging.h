#include "twoPhaseAveragingMacros.h"

double tAvStart = 0.0;
double tAvPrev = 0.0;

// Function which resets the averages by setting the start time and previous
// time to the current time

void resetAverages()
{
    tAvStart = t;
    tAvPrev = t;
}

// Divergence of a vector

void divergence(scalar d, vector u)
{
    trash({d});

    foreach()
    {
        d[] = 0.0;

        foreach_dimension()
        {
            d[] +=
                fm.x[1]*face_value(u.x,1)
              - fm.x[]*face_value(u.x,0);
        }

        d[] /= (Delta*cm[] + SEPS);
    }
}

// Gradient of a scalar

void grad(vector ga, scalar a)
{
    trash({ga});

    foreach()
    {
        foreach_dimension()
        {
            ga.x[] = center_gradient(a);
        }
    }
}

// Gradient of a scalar, with zero-gradient at the embedded wall

void gradzg(vector ga, scalar a)
{
    trash({ga});

    #if EMBED
    const double s = 1e-2;
    #endif

    foreach()
    {
        foreach_dimension()
        {
            #if EMBED

            if (fs.x[] > s && fs.x[1] > s)
            {
                ga.x[] = (a[1] - a[-1])/(2.0*Delta);
            }
            else if (fs.x[] > s)
            {
                ga.x[] = (a[] - a[-1])/(2.0*Delta);
            }
            else if (fs.x[1] > s)
            {
                ga.x[] = (a[1] - a[])/(2.0*Delta);
            }
            else
            {
                ga.x[] = 0.0;
            }

            #else

            ga.x[] = (a[1] - a[-1])/(2.0*Delta);

            #endif
        }
    }
}

// Laplacian of a scalar. Note that nu should contain the face metric, and in
// case of embedding the face aperture too. See
// https://groups.google.com/g/basilisk-fr/c/8qzXXaugIkY/m/2nkfsTV7CgAJ for the
// nowarning argument to the foreach call

void lapl(scalar la, face vector nu, scalar a)
{
    face vector fga[];

    foreach_face()
    {
        fga.x[] = nu.x[]*face_gradient_x(a,0);
    }

    foreach(nowarning)
    {
        la[] = 0.0;

        foreach_dimension()
        {
            la[] += (fga.x[1]-fga.x[])/Delta;
        }

        #if EMBED

        double c, d = embed_flux(point, a, fs, &c);

        if (c || d)
        {
            double nua = 0.0, fa = 0.0;

            foreach_dimension()
            {
                nua += nu.x[] + nu.x[1];
                fa  += fs.x[] + fs.x[1];
            }

            nua /= fa;

            la[] -= nua*(c+d*a[]);
        }

        #endif

        la[] /= (cm[] + SEPS);
    }
}

// Function to reconstruct vector from face vector

void reconstruct(vector a, face vector fa)
{
    trash({a});

    foreach()
    {
        foreach_dimension()
        {
            a.x[] =
                (fm.x[1]*fa.x[1] + fm.x[]*fa.x[])
              / (fm.x[1] + fm.x[] + SEPS);
        }
    }
}

// Interpolate data to a cylindrical mesh and average in circumferential
// direction. Write to the surf array.

void surfAverage
(
    double * surfs,
    scalar * fields,
    const double * aArr,
    const double * rArr,
    const double * tArr
)
{
    const int Nf = list_len(fields);

    // Find count

    int count[Na*Nr];

    for (int i = 0; i < Na; i++)
    for (int j = 0; j < Nr; j++)
    {
        count[j*Na+i] = 0;
    }

    // Initialize surfs to zero

    for (int l = 0; l < Nf; l++)
    for (int j = 0; j < Nr; j++)
    for (int k = 0; k < Na; k++)
    {
        surfs[l*Na*Nr+j*Na+k] = 0.0;
    }

    // Correct boundary conditions. Needed because the interpolation will not
    // automatically trigger it.

    boundary(fields);

    // Interpolate to structured grid and average in theta-direction

    for (int i = 0; i < Na; i++)
    {
        const double X = aArr[i];

        for (int j = 0; j < Nr; j++)
        {
            for (int k = 0; k < Nt; k++)
            {
                const double r = rArr[j];
                const double theta = tArr[k];

                const double Y = r*sin(theta);
                const double Z = r*cos(theta);

                Point point = locate(X, Y, Z);

                if (point.level >= 0)
                {
                    int l = 0;

                    count[j*Na+i]++;

                    for (scalar field in fields)
                    {
                        const double value =
                            interpolate_linear(point, (struct _interpolate){field, X, Y, Z});

                        surfs[l*Na*Nr+j*Na+i] += value/Nt;

                        l++;
                    }
                }
            }
        }
    }

    #if _MPI

    if (mpi_rank == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, surfs, Nf*Nr*Na, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, Nr*Na, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Reduce(surfs, surfs, Nf*Nr*Na, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, count, Nr*Na, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    #endif

    if (mpi_rank == 0)
        for (int i = 0; i < Na; i++)
            for (int j = 0; j < Nr; j++)
                if (count[j*Na+i] != Nt)
                    fprintf(stdout, "Could not locate all points for a = %g, r = %g %d %d\n", aArr[i], rArr[j], count[j*Na+i], Nt);

}

// Interpolate data to a Cartesian mesh and average in span- and streamwise
// direction. Write to the line array.

void lineAverage
(
    double * lines,
    scalar * fields,
    const double * xArr,
    const double * yArr,
    const double * zArr
)
{
    const int Nf = list_len(fields);

    // Find count

    int count[Nr];

    for (int j = 0; j < Nr; j++)
        count[j] = 0;

    // Initialize lines to zero

    for (int l = 0; l < Nf; l++)
    for (int j = 0; j < Nr; j++)
    {
        lines[l*Nr+j] = 0.0;
    }

    // Correct boundary conditions. Needed because the interpolation will not
    // automatically trigger it.

    boundary(fields);

    // Interpolate to Cartesian grid and average in x- and z-direction

    for (int i = 0; i < Na; i++)
    {
        const double X = xArr[i];

        for (int k = 0; k < Nt; k++)
        {
            const double Z = zArr[k];

            for (int j = 0; j < Nr; j++)
            {
                const double Y = yArr[j];

                Point point = locate(X, Y, Z);

                if (point.level >= 0)
                {
                    int l = 0;

                    count[j]++;

                    for (scalar field in fields)
                    {
                        const double value =
                            interpolate_linear(point, (struct _interpolate){field, X, Y, Z});

                        lines[l*Nr+j] += value/(Na*Nt);

                        l++;
                    }
                }
            }
        }
    }

    #if _MPI

    if (mpi_rank == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, lines, Nf*Nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, Nr, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Reduce(lines, lines, Nf*Nr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, count, Nr, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    #endif

    if (mpi_rank == 0)
        for (int j = 0; j < Nr; j++)
            if (count[j] != Na*Nt)
                fprintf(stdout, "Could not locate all points for y = %g\n", yArr[j]);
}

// Function to reconstruct the acceleration from a potential. Code replicates
// what the acceleration event in iforce.h does.

void reconstructAcceleration
(
    face vector a,
    scalar phi,
    scalar f
)
{
    foreach()
        f[] = clamp(f[], 0.0, 1.0);

    #if TREE
    f.prolongation = p.prolongation;
    f.dirty = true;
    #endif

    foreach_face()
    {
        if (f[] != f[-1])
        {
            double phif =

                (phi[] < nodata && phi[-1] < nodata)
              ? (phi[] + phi[-1])/2.0
              : phi[] < nodata ? phi[]
              : phi[-1] < nodata ? phi[-1]
              : 0.0;

            a.x[] += alpha.x[]/(fm.x[]+SEPS)*phif*(f[]-f[-1])/Delta;
        }
    }

    #if TREE
    f.prolongation = fraction_refine;
    f.dirty = true;
    #endif
}