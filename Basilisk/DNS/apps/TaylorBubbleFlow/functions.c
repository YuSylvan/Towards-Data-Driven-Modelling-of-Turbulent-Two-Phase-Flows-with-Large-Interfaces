#include <time.h>

// The recycle function assumes that the mesh remains unmodified throughout the
// simulation. Therefore, the points array can be cached

coord * points = NULL;
int numberOfPoints = 0;

#if _MPI
int * procOffsets = NULL;
#endif

void generateRecyclePoints()
{
    // Get boundary face counts

    int localNumberOfPoints = 0;

    foreach_boundary(right, reduction(+:numberOfPoints))
    {
        double r = 2.0*sqrt(sq(y)+sq(z));

        if (r <= 1.0+Delta)
        {
            localNumberOfPoints++;
            numberOfPoints++;
        }
    }

    // Store face coordinates

    double yl[localNumberOfPoints];

    #if dimension > 2
    double zl[localNumberOfPoints];
    #endif

    int k = 0;

    foreach_boundary(right)
    {
        double r = 2.0*sqrt(sq(y)+sq(z));

        if (r <= 1.0+Delta)
        {
            yl[k] = y;

            #if dimension > 2
            zl[k] = z;
            #endif

            k++;
        }
    }

    // If MPI, collect all face coordinates and create an array of points from
    // those coordinates. Otherwise, just create the array of points.

    points = malloc(sizeof(coord)*numberOfPoints);

    #if _MPI

        // Share face counts and compute processor offsets

        int counts[npe()];

        procOffsets = malloc(sizeof(int)*npe());

        MPI_Allgather
        (
            &localNumberOfPoints,
            1,
            MPI_INT,
            counts,
            1,
            MPI_INT,
            MPI_COMM_WORLD
        );

        procOffsets[0] = 0;

        for (int i = 1; i < npe(); i++)
        {
            procOffsets[i] = procOffsets[i-1]+counts[i-1];
        }

        // Collect all coordinates (MPI_Allgatherv doesn't work well for some
        // reason)

        double yg[numberOfPoints];

        MPI_Gatherv
        (
            yl,
            localNumberOfPoints,
            MPI_DOUBLE,
            yg,
            counts,
            procOffsets,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        MPI_Bcast
        (
            yg,
            numberOfPoints,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        #if dimension > 2

        double zg[numberOfPoints];

        MPI_Gatherv
        (
            zl,
            localNumberOfPoints,
            MPI_DOUBLE,
            zg,
            counts,
            procOffsets,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        MPI_Bcast
        (
            zg,
            numberOfPoints,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        #endif

        // Create points array from global coordinates

        for (int i = 0; i < numberOfPoints; i++)
        {
            points[i].x = X0+L0-SAMPLEDIST;
            points[i].y = yg[i];

            #if dimension > 2
            points[i].z = zg[i];
            #endif
        }

    #else

        // Create points array from local coordinates

        for (int i = 0; i < numberOfPoints; i++)
        {
            points[i].x = X0+L0-SAMPLEDIST;
            points[i].y = yl[i];

            #if dimension > 2
            points[i].z = zl[i];
            #endif
        }

    #endif
}

void recycle()
{
    // Sample the velocity field

    scalar * list = (scalar*){u};

    // Explicit boundary update needed for interpolation

    boundary(list);

    const int len = list_len(list);

    double * idata = malloc(sizeof(double)*len*numberOfPoints);

    interpolate_array(list, points, numberOfPoints, idata, false);

    double A = 0.0;
    double U = 0.0;

    // Broadcast sampled data back to all processors

    #if _MPI
    MPI_Bcast(idata, len*numberOfPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    // Set the sampled velocity field on the faces and return the bulk velocity

    int k = 0;

    foreach_boundary(right, reduction(+:A) reduction(+:U))
    {
        double r = 2.0*sqrt(sq(y)+sq(z));

        if (r <= 1.0+Delta)
        {
            #if _MPI
            const int i = (procOffsets[mpi_rank] + k)*len;
            #else
            const int i = k;
            #endif

            v.x[] = idata[i];
            v.y[] = idata[i+1];
            v.z[] = idata[i+2];

            const double Ai = sq(Delta)*fs.x[];

            A += Ai;
            U += Ai*v.x[];

            k++;
        }
    }

    U /= A;

    // Adjust sampled velocity field to have the right bulk velocity

    foreach_boundary(right)
    {
        v.x[] = v.x[] - U + Ui;
    }

    free(idata);
}

void updateApertures()
{
    vertex scalar phi[];

    foreach_vertex()
    {
        phi[] = sq(0.5*1.000001) - sq(y) - sq(z);
    }

    fractions(phi, cs, fs);
}

double liquidVelocity(double xp)
{
    double A = 0.0;
    double flowRate = 0.0;

    foreach_face
    (
        x,
        reduction(+:A)
        reduction(+:flowRate)
    )
    {
        if (x > xp-Delta/2.0 && x <= xp+Delta/2.0)
        {
            A += sq(Delta)*fs.x[];

            flowRate += sq(Delta)*face_value(f,0)*uf.x[];
        }
    }

    return flowRate/A;
}

double bubbleVelocity(double xp)
{
    double A = 0.0;
    double flowRate = 0.0;

    foreach_face
    (
        x,
        reduction(+:A)
        reduction(+:flowRate)
    )
    {
        if (x > xp-Delta/2.0 && x <= xp+Delta/2.0)
        {
            A += sq(Delta)*fs.x[];

            flowRate += sq(Delta)*(1.0-face_value(f,0))*uf.x[];
        }
    }

    return flowRate/A;
}

double totalVoidFraction()
{
    double vol = 0.0;
    double totalf = 0.0;

    foreach
    (
        reduction(+:vol)
        reduction(+:totalf)
    )
    {
        vol += dv()*cs[];
        totalf += (1.0-f[])*dv()*cs[];
    }

    return totalf/vol;
}

void initLiquid()
{
    foreach()
    {
        f[] = 1.0;
    }
}

void initPlugVelocity(const double U)
{
    foreach()
    {
        u.x[] = cs[]*U;
        u.y[] = 0.0;
        u.z[] = 0.0;
    }
}

void initBubble()
{
    foreach()
    {
        if
        (
            (x > XB0-LB0/2.0 && x <= XB0+LB0/2.0 && sq(y) + sq(z) < sq(DB0/2.0))
         || (sq(x-(XB0+LB0/2.0)) + sq(y) + sq(z) < sq(DB0/2.0))
        )
        {
            f[] = 0.0;
        }
        else
        {
            f[] = 1.0;
        }
    }

    // Place the bubble in a zero velocity region that extends a quarter pipe
    // diameter below and above the bubble.

    foreach()
    {
        if (x > XB0-LB0/2.0-0.25 && x <= XB0+LB0/2.0+DB0/2.0+0.25)
        {
            foreach_dimension()
            {
                u.x[] = 0.0;
            }
        }
    }
}

void updateFlux()
{
    trash ({uf});

    foreach_face()
    {
        uf.x[] = fs.x[]*face_value(u.x, 0);
    }
}

void scaleAndShiftVelocity(double U)
{
    // It is assumed here that a single phase turbulent flow field is provided
    // as generated by the pipeFlow app. Rescale the velocity by the given
    // factor and subtract offset to give zero mean.

    foreach()
    {
        u.x[] = u.x[]/U - 1.0;
        u.y[] = u.y[]/U;
        u.z[] = u.z[]/U;
    };
}

double bubbleVolume(coord * p)
{
    scalar m[];

    foreach()
    {
        m[] = f[] < 0.999;
    }

    int n = tag(m);

    double v[n];
    coord b[n];

    for (int j = 0; j < n; j++)
    {
        v[j] = b[j].x = b[j].y = b[j].z = 0.0;
    }

    foreach_leaf()
    {
        if (m[] > 0)
        {
            int j = m[] - 1;

            v[j] += dv()*(1.0-f[]);
            coord p = {x,y,z};

            foreach_dimension()
            {
                b[j].x += dv()*(1.0-f[])*p.x;
            }
        }
    }

    #if _MPI
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    // It is assumed that the largest volume is the Taylor bubble

    int maxj = 0;

    for (int j = 0; j < n; j++)
    {
        if (v[j] >= v[maxj])
        {
            maxj = j;
        }
    }

    p->x = b[maxj].x/v[maxj];
    p->y = b[maxj].y/v[maxj];
    p->z = b[maxj].z/v[maxj];

    return v[maxj];
}
