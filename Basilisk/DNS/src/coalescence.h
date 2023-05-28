// This code models coalescence between particles. Important: particles are
// defined by unity volume fraction f. So the dispersed phase (either liquid or
// gas) should be associated with f[] = 1.0. So switching density/viscosity may
// be needed.

// Two parameters must be set:
//
//     FILM_THICKNESS: the distance, scaled by cell size, between which the
//     centroids of two interface cell fractions must lie to be considered as
//     enclosing a film. Defaults to 1.
//
//     CRITICAL_TIME: the dimensional time that must be surpassed in order to
//     allow coalescence. This quantity is currently also a constant and should
//     also be modeled in the future

#ifndef FILM_THICKNESS
#define FILM_THICKNESS 1.0
#endif

#include "tag.h"

#define BGHOSTS 2

// Field to collect color numbers in

scalar colors[];

// For each interface a tag field

scalar * tags = NULL;

// Global number of particles

int Np = 0;
int NpOld = 0;

// For each particle the interface number

int * interfaceNums;

// Collision pair data

typedef struct
{
    int p1;
    int p2;
    double t;
    double We;
    double xi;
    double A1;
    double A2;
} collisionPair;

collisionPair emptyPair = {-1, -1, 0, 0, 0, 0, 0};

collisionPair * pairs;

int Npairs = 0;

// Variable to disable surface tension

#ifndef NOTENSION
#define NOTENSION 0
#endif

// Function to add an interface to the list of interfaces

void addInterface()
{
    // Create new field on te heap

    scalar fnew = new scalar;

    foreach()
    {
        fnew[] = 0.0;
    }

    #if !NOTENSION
    fnew.sigma = f.sigma;
    #endif

    #if TREE
    fnew.restriction = restriction_volume_average;
    fnew.refine = fraction_refine;
    fnew.prolongation = fraction_refine;
    fnew.dirty = true;
    #endif

    // Update the name of the field

    char * name = malloc(sizeof(char)*100);
    sprintf(name, "f%d", list_len(interfaces));

    free(fnew.name);
    fnew.name = name;

    // Append to interfaces

    interfaces = list_append(interfaces,fnew);

    fprintf
    (
        stdout,
        "Coalescence: Added %s to interfaces. Number of interfaces = %d\n",
        name,
        list_len(interfaces)
    );
}

// Function to add a tag field to the list of tag fields

void addTag()
{
    scalar g = new scalar;

    // Tag numbers should not be averaged/interpolated in some way. They should
    // be injected, which preserves round numbers.

    #if TREE
    g.restriction = restriction_tag;
    g.prolongation = refine_injection;
    g.refine = refine_injection;
    g.dirty = true;
    #endif

    foreach()
    {
        g[] = 0.0;
    }

    // Update the name of the field

    char * name = malloc(sizeof(char)*100);
    sprintf(name, "g%d", list_len(tags));

    free(g.name);
    g.name = name;

    // Append to the list of tags

    tags = list_append(tags,g);

    fprintf
    (
        stdout,
        "Coalescence: Added %s to tags. Number of tags = %d\n",
        name,
        list_len(tags)
    );
}

// Function to store the collision pairs

void storeCollisionPairs
(
    const char fileName[],
    const int dirIndex
)
{
    // Only the master process needs to do this

    #if _MPI

    if (mpi_rank > 0)
        return;

    #endif

    char dirName[20];

    sprintf(dirName, "%06d", dirIndex);

    DIR * dir;
    dir = NULL;

    if (access(dirName, F_OK) == -1 || (dir = opendir(dirName)) == NULL)
    {
        mkdir(dirName, 0755);
    }

    if (dir)
    {
        closedir(dir);
    }

    FILE * file;

    char filePath[80];

    sprintf(filePath, "%s/%s", dirName, fileName);

    file = fopen(filePath, "w");

    // Write the number of particles and pairs

    fprintf(file, "%d\n", Np);
    fprintf(file, "%d\n", Npairs);

    // Write the collision pairs on separate lines

    for (int i = 0; i < Npairs; i++)
        fprintf
        (
            file,
            "%d %d %g %g %g %g %g\n",
            pairs[i].p1,
            pairs[i].p2,
            pairs[i].t,
            pairs[i].We,
            pairs[i].xi,
            pairs[i].A1,
            pairs[i].A2
        );

    fclose(file);
}

// Function to restore tags, interfaces and collision pairs. Fields are restored
// by adding all "f[0-9]+" and "f[0-9]+" fields that are registered to cleared
// interface and tag lists. Needed when resuming a simulation.

void restoreCoalescence
(
    const char fileName[]
)
{
    interfaces = NULL;
    tags = NULL;

    char name[100];

    // Restore interfaces and tags by searching through all registered fields.
    // It is important that the restore command is called with the "list = all"
    // option, so that fields are restored even though they are not present yet.

    for (scalar s in all)
    {
        for (int i = 0; i < 100; i++)
        {
            sprintf(name, "f%d", i);

            if (strcmp(s.name,name) == 0)
            {
                interfaces = list_append(interfaces, s);

                #if !NOTENSION
                s.sigma = f.sigma;
                #endif

                #if TREE
                s.restriction = restriction_volume_average;
                s.refine = fraction_refine;
                s.prolongation = fraction_refine;
                s.dirty = true;
                #endif

                fprintf
                (
                    stdout,
                    "Coalescence: Restored %s to interfaces\n",
                    name
                );
            }

            sprintf(name, "g%d", i);

            if (strcmp(s.name,name) == 0)
            {
                tags = list_append(tags, s);

                #if TREE
                s.restriction = restriction_tag;
                s.prolongation = refine_injection;
                s.refine = refine_injection;
                s.dirty = true;
                #endif

                fprintf
                (
                    stdout,
                    "Coalescence: Restored %s to tags\n",
                    name
                );
            }
        }
    }

    fprintf
    (
        stdout,
        "Coalescence: Reset interfaces and tags. "
        "Number of interfaces = %d, number of tags = %d\n",
        list_len(interfaces),
        list_len(tags)
    );

    // Restore collision pairs

    if (access(fileName, F_OK) != -1)
    {
        FILE* file;

        file = fopen(fileName, "r");

        fscanf(file, "%d\n", &Np);
        fscanf(file, "%d\n", &Npairs);

        pairs = malloc(Np*sizeof(collisionPair));

        char buffer[1024];
        char * ptr;
        int i = 0;

        while ((ptr = fgets(buffer, sizeof(buffer), file)) != NULL)
        {
            sscanf
            (
                buffer,
                "%d %d %lf %lf %lf %lf %lf\n",
                &pairs[i].p1,
                &pairs[i].p2,
                &pairs[i].t,
                &pairs[i].We,
                &pairs[i].xi,
                &pairs[i].A1,
                &pairs[i].A2
            );

            i++;
        }

        fclose(file);
    }
    else
    {
        fprintf
        (
            stdout,
            "Coalescence: Could not find coalescence restore file %s. " "Aborting\n",
            fileName
        );
        exit(1);
    }
}

// Function to collect all color numbers into a single color field. Useful for
// visualization.

void collectColors()
{
    foreach()
    {
        colors[] = 0.0;
    }

    int interfaceNum = 0;

    for(scalar fi in interfaces)
    {
        foreach()
        {
            colors[] += round(fi[])*(interfaceNum+1);
        }

        interfaceNum++;
    }
}

// Function to collect all interfaces into f

void collectInterfaces()
{
    foreach()
    {
        f[] = 0.0;
    }

    for (scalar fi in interfaces)
    {
        foreach()
        {
            f[] = min(f[]+fi[],1.0);
        }
    }
}

// Function that sets the particle interface number array given the number of
// particles per interface.

void setInterfaceNums(int * Nps)
{
    // Reset the length of the interface number array based on the total number
    // of particles

    free(interfaceNums);
    interfaceNums = malloc(Np*sizeof(int));

    // Fill the interface array

    int k = 0;

    for (int i = 0; i < list_len(interfaces); i++)
        for (int j = 0; j < Nps[i]; j++)
            interfaceNums[k++] = i;
}

// Function that initializes and updates the connectivity matrix, by finding
// close particle pairs.

void setConnectivityMatrix(bool *** CPtr)
{
    // Allocate the 2d connectivity matrix contiguously, so that the matrix can
    // be updated in one MPI call.

    *CPtr = malloc(Np*sizeof(bool*));

    bool * data = malloc(Np*Np*sizeof(bool));

    for (int i = 0; i < Np; i++)
    {
        (*CPtr)[i] = data + Np*i;
    }

    bool ** C = *CPtr;

    for (int i = 0; i < Np; i++)
    for (int j = 0; j < Np; j++)
    {
        C[i][j] = false;
    }

    // Create the shadow field, collecting all particle numbers into one field,
    // from which connectivity can be determined. On trees, assure that
    // interpolations are preserving round numbers, by injecting.

    scalar shadow[];

    #if TREE
    shadow.restriction = restriction_tag;
    shadow.prolongation = refine_injection;
    shadow.refine = refine_injection;
    shadow.dirty = true;
    #endif

    foreach()
    {
        shadow[] = 0.0;

        for (scalar g in tags)
        {
            if (round(shadow[]) < 1.0 && round(g[]) > 0.0)
            {
                shadow[] = g[];
            }
        }
    }

    // Fill the connectivity matrix

    #if dimension == 2
    int Nh = 25;
    #else
    int Nh = 125;
    #endif

    int tagNums[Nh];

    for (int i = 0; i < Nh; i++)
        tagNums[i] = 0;

    foreach()
    {
        // In a box of size 5 grid cells around the current grid cell, find all
        // distinct tag numbers and store in a list

        int i = 0;

        foreach_neighbor()
        {
            if (shadow[] > 0.0)
            {
                int tagNum = (int) round(shadow[]);

                bool found = false;

                for (int j = 0; j < i; j++)
                if (tagNums[j] == tagNum)
                {
                    found = true;
                    break;
                }

                if (!found)
                    tagNums[i++] = tagNum;
            }
        }

        // For all combinations of found distinct tag numbers, set the
        // connectivity matrix

        for (int j = 0; j < i-1; j++)
        for (int k = j+1; k < i; k++)
        {
            C[tagNums[j]-1][tagNums[k]-1] = true;
            C[tagNums[k]-1][tagNums[j]-1] = true;
        }
    }

    #if _MPI

    // Communicate connectivity by a collective 'or' operation. Act on the
    // contiguous data, by treating the first element of C as an Np^2
    // one-dimensional array. This reduces the number of MPI calls to just one.

    MPI_Allreduce(MPI_IN_PLACE, *C, Np*Np, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    #endif
}

// Function to move a particle (identified by its global number) to an interface
// which has no neighboring particles. If no such interface exists, create a new
// one.

void moveParticle(bool ** C, int particleNum)
{
    // Particle interface number

    const int interfaceNum = interfaceNums[particleNum];

    // Create ineligible interface array. This array identifies which interfaces
    // are ineligible for hosting the particle.

    bool ineligible[list_len(interfaces)];

    for (int i = 0; i < list_len(interfaces); i++)
    {
        ineligible[i] = false;
    }

    // The current interface number is ineligible, by definition

    ineligible[interfaceNum] = true;

    // Fields that have particles that are connected to the current one are
    // ineligible too

    for (int i = 0; i < Np; i++)
    {
        int interfaceNum2 = interfaceNums[i];

        if (interfaceNum != interfaceNum2 && C[i][particleNum])
        {
            ineligible[interfaceNum2] = true;
        }
    }

    // Communicate the ineligible array by a collective 'or' operation

    #if _MPI
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        ineligible,
        list_len(interfaces),
        MPI_C_BOOL,
        MPI_LOR,
        MPI_COMM_WORLD
    );
    #endif

    // Find eligible interface

    int targetNum = -1;

    for (int i = 0; i < list_len(interfaces); i++)
    {
        if (!ineligible[i])
        {
            targetNum = i;
        }
    }

    // Create and append new interface if there is no eligible one. The target
    // interface becomes this newly created one.

    if (targetNum < 0)
    {
        addInterface();
        addTag();

        targetNum = list_len(interfaces)-1;
    }

    // Move the particle from the source interface to the target interface

    const int sourceNum = interfaceNums[particleNum];

    scalar source = interfaces[sourceNum];
    scalar target = interfaces[targetNum];

    scalar gSource = tags[sourceNum];
    scalar gTarget = tags[targetNum];

    foreach()
    {
        if (gSource[] > 0.0 && round(gSource[]-1.0) == particleNum)
        {
            gTarget[] = particleNum + 1.0;
            gSource[] = 0.0;

            target[] = min(source[]+target[],1.0);
            source[] = 0.0;
        }
    }

    // Update interface number

    interfaceNums[particleNum] = targetNum;

    // Update connectivity matrix. By definition, the particle is moved in such
    // a way that it is no longer connected to any other particle, so the
    // particleNum's row and column in the connectivity matrix should be false.

    for (int i = 0; i < Np; i++)
    {
        C[particleNum][i] = false;
        C[i][particleNum] = false;
    }

    fprintf
    (
        stdout,
        "Coalescence: Moved particle %d from interface %d to interface %d\n",
        particleNum,
        sourceNum,
        targetNum
    );
}

// Function that finds connected pairs and then moves them

void moveConnectedPairs(bool ** C)
{
    // Sweep only the lower/upper triangle of the connectivity matrix to avoid
    // redundancy

    for (int i = 0; i < Np-1; i++)
    {
        for (int j = i+1; j < Np; j++)
        {
            // Relocate second particle if both particles are in the same
            // interface and are connected

            if (C[i][j] && interfaceNums[i] == interfaceNums[j])
            {
                fprintf
                (
                    stdout,
                    "Coalescence: Found connected particle pair (%d,%d) in interface %d\n",
                    i,
                    j,
                    interfaceNums[i]
                );

                moveParticle(C,j);
            }
        }
    }
}

// Function that prevents coalescence by separating close particle pairs

void preventCoalescence()
{
    // Create the connectivity matrix

    bool ** C = NULL;
    setConnectivityMatrix(&C);

    // Move close particle pairs to different interfaces, to prevent
    // numerical coalescence

    moveConnectedPairs(C);

    // Free up memory

    free(C[0]);
    free(C);
}

// Function to calculate particle volumes. In 2d this returns particle area
// because that's what dv() does.

void calcParticleVolumes(double ** VPtr)
{
    *VPtr = malloc(Np*sizeof(double));

    double * V = *VPtr;

    for (int i = 0; i < Np; i++)
        V[i] = 0.0;

    for (int i = 0; i < list_len(interfaces); i++)
    {
        scalar f = interfaces[i];
        scalar g = tags[i];

        foreach()
        {
            if (g[] > 0)
            {
                int particleNum = (int) round(g[]-1.0);
                V[particleNum] += f[]*dv();
            }
        }
    }

    #if _MPI
    MPI_Allreduce(MPI_IN_PLACE, V, Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
}

// Function to calculate particle velocity vectors

void calcParticleVelocities(coord ** UPtr, double * V)
{
    *UPtr = malloc(Np*sizeof(coord));

    coord * U = *UPtr;

    for (int i = 0; i < Np; i++)
    {
        U[i].x = 0;
        U[i].y = 0;
        U[i].z = 0;
    }

    // Volume-integrate velocity field for each particle

    for (int i = 0; i < list_len(interfaces); i++)
    {
        scalar f = interfaces[i];
        scalar g = tags[i];

        foreach()
        {
            if (g[] > 0)
            {
                int particleNum = (int) round(g[]-1.0);

                U[particleNum].x += u.x[]*f[]*dv();
                U[particleNum].y += u.y[]*f[]*dv();
                U[particleNum].z += u.z[]*f[]*dv();
            }
        }
    }

    #if _MPI
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        U,
        Np*sizeof(coord)/sizeof(double),
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    #endif

    // Scale by particle volume

    for (int i = 0; i < Np; i++)
    {
        U[i].x /= V[i];
        U[i].y /= V[i];
        U[i].z /= V[i];
    }
}

// Function to compute the total interface between two particles from a cloud of
// fragments. The function returns a fragment with a barycenter and the total
// area of the interface.

typedef struct { coord c; double a; } fragment;

fragment pairInterface(fragment * Fs, int N)
{
    fragment F = {{0.0, 0.0, 0.0}, 0.0};

    for (int i = 0; i < N; i++)
    {
        F.a += Fs[i].a;

        F.c.x += Fs[i].c.x * Fs[i].a;
        F.c.y += Fs[i].c.y * Fs[i].a;
        F.c.z += Fs[i].c.z * Fs[i].a;
    }

    F.c.x /= F.a;
    F.c.y /= F.a;
    F.c.z /= F.a;

    // To debug fragments:

    // fprintf
    // (
    //     stdout,
    //     "N = %d, c = (%g, %g, %g), a = %g:\n",
    //     N,
    //     F.c.x,
    //     F.c.y,
    //     F.c.z,
    //     F.a
    // );

    // for (int i = 0; i < N; i++)
    // {
    //     fprintf
    //     (
    //         stdout,
    //         "    i = %d, c = (%g, %g, %g), a = %g\n",
    //         i,
    //         Fs[i].c.x,
    //         Fs[i].c.y,
    //         Fs[i].c.z,
    //         Fs[i].a
    //     );
    // }

    return F;
}

// Function to calculate the interface areas matrix

void calcAreaMatrix(double *** APtr, double * V)
{
    // Allocate the 2d area matrix contiguously, so that the matrix can
    // be updated in one MPI call.

    *APtr = malloc(Np*sizeof(double*));

    double * data = malloc(Np*Np*sizeof(double));

    for (int i = 0; i < Np; i++)
    {
        (*APtr)[i] = data + Np*i;
    }

    double ** A = *APtr;

    for (int i = 0; i < Np; i++)
    for (int j = 0; j < Np; j++)
    {
        A[i][j] = 0.0;
    }

    // For all interface combinations, find close fragments contained by two
    // distinct fields. Store information on each interface fragment in a list.
    // A fragment contains position and area.

    // Prepare fragment areas and centers for each field

    scalar * areas = NULL;
    scalar * centers = NULL;

    for (int i = 0; i < list_len(interfaces); i++)
    {
        scalar area = new scalar;
        vector center = new vector;

        scalar f = interfaces[i];
        scalar g = tags[i];

        foreach()
        if (f[] > 0.01 && f[] < 0.99 && g[] > 0.0)
        {
            coord n = interface_normal(point, f);
            double alpha = plane_alpha(f[], n);

            coord c;

            #if dimension == 2
            area[] = line_length_center(n, alpha, &c)*Delta;
            #else
            area[] = plane_area_center(n, alpha, &c)*sq(Delta);
            #endif

            center.x[] = x + c.x*Delta;
            center.y[] = y + c.y*Delta;
            center.z[] = z + c.z*Delta;
        }
        else
        {
            area[] = 0.0;

            center.x[] = 0.0;
            center.y[] = 0.0;
            center.z[] = 0.0;
        }

        areas = list_append(areas,area);
        centers = list_concat(centers, (scalar*){center});
    }

    // Create fragment lists

    fragment * fragmentListPtrs[Np*Np];

    int Ns[Np*Np];

    for (int i = 0; i < Np*Np; i++)
    {
        fragmentListPtrs[i] = NULL;
        Ns[i] = 0;
    }

    // Fill fragment lists

    const int incr = 16;

    for (int i = 0; i < list_len(interfaces); i++)
    for (int j = 0; j < list_len(interfaces); j++)
    if (i != j)
    {
        scalar gi = tags[i];
        scalar gj = tags[j];

        scalar ai = areas[i];
        scalar aj = areas[j];

        #if dimension == 2
        vector ci = {centers[2*i], centers[2*i+1]};
        vector cj = {centers[2*j], centers[2*j+1]};
        #else
        vector ci = {centers[3*i], centers[3*i+1], centers[3*i+2]};
        vector cj = {centers[3*j], centers[3*j+1], centers[3*j+2]};
        #endif

        foreach()
        if (ai[] > 0.0)
        {
            double area = ai[];
            coord center = {ci.x[], ci.y[], ci.z[]};

            int a = (int) round(gi[]-1.0);

            // Search in a stencil of width 3 to see if the other interface also
            // has a fragment close by

            bool found = false;

            foreach_neighbor(1)
            if (!found && aj[] > 0.0)
            {
                // If the interface fragment centroids are within the
                // non-dimensional film thickness threshold, add.

                if
                (
                    sq(center.x - cj.x[])
                  + sq(center.y - cj.y[])
                  + sq(center.z - cj.z[])
                  < sq(FILM_THICKNESS*Delta)
                )
                {
                    // Particle number in interface j

                    int b = (int) round(gj[]-1.0);

                    // Index of pair (a,b) in the matrix A

                    int c = a*Np + b;

                    // Append fragment to the list

                    if (Ns[c] == 0)
                    {
                        // Create buffer

                        fragmentListPtrs[c] =
                            malloc(incr*sizeof(fragment));
                    }
                    else if ((Ns[c] % incr) == 0)
                    {
                        // Increase buffer size

                        fragmentListPtrs[c] = realloc
                        (
                            fragmentListPtrs[c],
                            (Ns[c]+incr)*sizeof(fragment)
                        );
                    }

                    // Add fragment

                    fragment F = {center, area};

                    fragmentListPtrs[c][Ns[c]] = F;

                    // Once the interface i fragment was found to be close to an
                    // interface j fragment, we don't have to look further.

                    found = true;

                    Ns[c]++;
                }
            }
        }
    }

    // Free memory

    delete(areas);
    delete(centers);

    free(areas);
    free(centers);

    // Calculate the array of pair interfaces

    int Ni = 0, p = 0;

    #if _MPI

    int Nsa[Np*Np];

    MPI_Allreduce(Ns, Nsa, Np*Np, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < Np*Np; i++)
        if (Nsa[i] > 0)
            Ni++;

    fragment pairInterfaces[Ni];

    for (int i = 0; i < Np*Np; i++)
    {
        if (Nsa[i] > 0)
        {
            // For this pair, gather all processor's Ns values to master

            int Nsp[npe()];

            MPI_Gather(&Ns[i], 1, MPI_INT, Nsp, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Now gather all fragments (as doubles) to master

            int disp[npe()];
            int size[npe()];

            disp[0] = 0;

            for (int j = 1; j < npe(); j++)
                disp[j] = disp[j-1] + Nsp[j-1]*4;

            for (int j = 0; j < npe(); j++)
                size[j] = Nsp[j]*4;

            fragment * Fs = malloc(Nsa[i]*sizeof(fragment));

            MPI_Gatherv
            (
                fragmentListPtrs[i],
                Ns[i]*4,
                MPI_DOUBLE,
                Fs,
                size,
                disp,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD
            );

            // Get the pair interface from this cloud of fragments

            if (mpi_rank == 0)
                pairInterfaces[p] = pairInterface(Fs, Nsa[i]);

            p++;
        }
    }

    // Broadcast all pair interfaces

    MPI_Bcast(pairInterfaces, Ni*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #else

    for (int i = 0; i < Np*Np; i++)
        if (Ns[i] > 0)
            Ni++;

    fragment pairInterfaces[Ni];

    for (int i = 0; i < Np*Np; i++)
    {
        if (Ns[i] > 0)
        {
            // Get the pair interface from this cloud of fragments

            pairInterfaces[p] = pairInterface(fragmentListPtrs[i], Ns[i]);

            p++;
        }
    }

    #endif

    // Calculate non-dimensional areas

    p = 0;

    for (int i = 0; i < Np*Np; i++)
    {
        #if _MPI
        if (Nsa[i] > 0)
        #else
        if (Ns[i] > 0)
        #endif
        {
            const int a = i/Np;
            const int b = i%Np;

            // In 2d V is an area while in 3d it's a volume. Estimate the
            // particle radius from a circle in 2d or sphere in 3d.

            #if dimension == 2
            double ra = sqrt(V[a]/pi);
            double rb = sqrt(V[b]/pi);
            #else
            double ra = cbrt(V[a]*3.0/pi/4.0);
            double rb = cbrt(V[b]*3.0/pi/4.0);
            #endif

            // In 2d the fragment's area is a length while in 3d it's an actual
            // area. Estimate the pair interface radius from a line in 2d or a
            // circle in 3d.

            #if dimension == 2
            double R = pairInterfaces[p].a/2.0;
            #else
            double R = sqrt(pairInterfaces[p].a/pi);
            #endif

            // Store the non-dimensional area

            A[a][b] = sq(R/ra);
            A[b][a] = sq(R/rb);

            // To debug interface areas:

            // fprintf
            // (
            //     stdout,
            //     "Interface pair (%d,%d), area = (%g,%g)\n",
            //     a,
            //     b,
            //     sq(R/ra),
            //     sq(R/rb)
            // );

            p++;
        }
    }

    // Free memory

    for (int i = 0; i < Np*Np; i++)
    {
        if (Ns[i] > 0)
        {
            free(fragmentListPtrs[i]);
        }
    }
}

// Function to merge a given pair of particles by moving the second particle
// into the first particle's interface

void mergePair(int p)
{
    collisionPair * pair = &pairs[p];

    const int sourceNum = interfaceNums[pair->p2];
    const int targetNum = interfaceNums[pair->p1];

    scalar source = interfaces[sourceNum];
    scalar target = interfaces[targetNum];

    scalar gSource = tags[sourceNum];
    scalar gTarget = tags[targetNum];

    // Move void

    foreach()
    {
        if (gSource[] > 0.0 && round(gSource[]-1.0) == pair->p2)
        {
            gTarget[] = pair->p1+1.0;
            gSource[] = 0.0;

            target[] = min(source[]+target[],1.0);
            source[] = 0.0;
        }
    }

    // All pairs that were in contact with p2 are now in contact with p1

    for (int i = 0; i < Npairs; i++)
    if (i != p)
    {
        if (pairs[i].p1 == pair->p2)
        {
            pairs[i].p1 = pair->p1;
        }
        else if (pairs[i].p2 == pair->p2)
        {
            pairs[i].p2 = pair->p1;
        }

        // TODO: Do something about We, xi, A1 and A2?
    }

    fprintf
    (
        stdout,
        "Coalescence: Merged particle %d into %d\n",
        pair->p2,
        pair->p1
    );
}

// Function that merges pairs based on a critical interface area

void mergePairs()
{
    for (int p = 0; p < Npairs; p++)
    {
        // TODO: Some sensible criterion here. It is based on a constant
        // TODO: critical value now. The collisionPair struct returned by
        // TODO: pairs[p] provides xi and We.

        if (pairs[p].t > CRITICAL_TIME)
        {
            mergePair(p);
        }
    }
}

// Function to transfer a tag field while keeping track of particle number
// conversion

void transferTag
(
    scalar gOld,
    scalar gNew,
    int * oldToNewNums
)
{
    // Get largest and smallest tag numbers in old and new tag field

    int minOld = 9999;
    int minNew = 9999;
    int maxOld = 0;
    int maxNew = 0;

    foreach
    (
        reduction(min:minOld)
        reduction(min:minNew)
        reduction(max:maxOld)
        reduction(max:maxNew)
    )
    {
        if (gOld[] > 0.0)
        {
            int tag = round(gOld[]);

            minOld = min(minOld,tag);
            maxOld = max(maxOld,tag);
        }

        if (gNew[] > 0.0)
        {
            int tag = round(gNew[]);

            minNew = min(minNew,tag);
            maxNew = max(maxNew,tag);
        }
    }

    // Number of tags in the old and new field

    int nOld = max(maxOld-minOld+1,0);
    int nNew = max(maxNew-minNew+1,0);

    if (nOld > 0 && nNew > 0)
    {
        // We have both old and new particle tags. Let's match them by
        // volumetric overlap.

        double crossVolume[nOld][nNew];
        double oldVolume[nOld];

        for (int i = 0; i < nOld; i++)
            for (int j = 0; j < nNew; j++)
                crossVolume[i][j] = 0.0;

        for (int i = 0; i < nOld; i++)
            oldVolume[i] = 0.0;

        foreach()
        {
            if (gOld[] > 0.0)
            {
                int i = round(gOld[]) - minOld;

                oldVolume[i] += dv();

                if (gNew[] > 0.0)
                {
                    int j = round(gNew[]) - minNew;

                    crossVolume[i][j] += dv();
                }
            }
        }

        #if _MPI

        MPI_Allreduce
        (
            MPI_IN_PLACE,
            crossVolume[0],
            nOld*nNew,
            MPI_DOUBLE,
            MPI_SUM,
            MPI_COMM_WORLD
        );

        MPI_Allreduce
        (
            MPI_IN_PLACE,
            oldVolume,
            nOld,
            MPI_DOUBLE,
            MPI_SUM,
            MPI_COMM_WORLD
        );

        #endif

        // Fill old-to-new tag array

        for (int i = 0; i < nOld; i++)
        {
            for (int j = 0; j < nNew; j++)
            {
                // If we have 60% volumetric overlap, consider the tags to
                // correspond to the same particle

                if (crossVolume[i][j]/max(oldVolume[i],1e-12) > 0.6)
                {
                    if (oldToNewNums[i+minOld-1] != -1)
                    {
                        fprintf
                        (
                            stdout,
                            "Coalescence: Old particle %d overlaps two "
                            "new particles (%d and %d)\n",
                            i+minOld-1,
                            oldToNewNums[i+minOld-1],
                            j+minNew-1
                        );
                    }

                    oldToNewNums[i+minOld-1] = j+minNew-1;
                }
            }
        }
    }

    // Transfer tags

    foreach()
    {
        gOld[] = gNew[];
    }
}

// Function to check the old-to-new particle number conversion array

void checkOldToNewNums
(
    int * oldToNewNums,
    int NpOld,
    int NpNew
)
{
    // Look for unmatched particles or doubly matched particles

    int match[Np];

    for(int i = 0; i < NpNew; i++)
        match[i] = -1;

    for(int i = 0; i < NpOld; i++)
    {
        int j = oldToNewNums[i];

        if (j == -1)
        {
            fprintf
            (
                stdout,
                "Coalescence: Particle %d could not be matched with a new "
                "particle tag. NpNew = %d and NpOld = %d\n",
                i,
                NpNew,
                NpOld
            );
        }
        else
        {
            if (match[j] != -1)
            {
                fprintf
                (
                    stdout,
                    "Coalescence: New particle %d overlaps two old particles "
                    "(%d and %d)\n",
                    j,
                    i,
                    match[j]
                );
            }
            else
            {
                match[j] = i;
            }
        }
    }
}

// Function to update collision pairs

void updateCollisionPairs(int * oldToNewNums)
{
    // Update particle numbers from the particle number conversion array. For
    // old particles that were not matched with a new one, the p1 or p2 value
    // will become -1.

    for (int i = 0; i < Npairs; i++)
    {
        pairs[i].p1 = oldToNewNums[pairs[i].p1];
        pairs[i].p2 = oldToNewNums[pairs[i].p2];
    }

    // Calculate particle volumes

    double * V = NULL;
    calcParticleVolumes(&V);

    // Calculate interface areas

    double ** A = NULL;
    calcAreaMatrix(&A, V);

    // Count new pairs

    int NpairsNew = 0;

    for (int i = 0; i < Np-1; i++)
        for (int j = i+1; j < Np; j++)
            if (A[i][j] > 0)
                NpairsNew++;

    // Calculate particle velocities, needed for the Weber number

    coord * U = NULL;
    calcParticleVelocities(&U, V);

    // Create and fill new pair array. Pairs that are no longer in contact
    // according to the area matrix are removed automatically, because they are
    // not copied into the new array.

    collisionPair * pairsNew = malloc(NpairsNew*sizeof(collisionPair));

    for (int i = 0; i < NpairsNew; i++)
        pairsNew[i] = emptyPair;

    int c = 0;

    for (int i = 0; i < Np-1; i++)
    {
        for (int j = i+1; j < Np; j++)
        {
            if (A[i][j] > 0)
            {
                // Look in old pairs array. Copy and update if found. If not,
                // create new pair.

                bool found = false;

                for(int p = 0; p < Npairs; p++)
                {
                    if
                    (
                        (pairs[p].p1 == i && pairs[p].p2 == j)
                     || (pairs[p].p1 == j && pairs[p].p2 == i)
                    )
                    {
                        pairsNew[c] = pairs[p];

                        pairsNew[c].t += dt;

                        pairsNew[c].A1 = A[i][j];
                        pairsNew[c].A2 = A[j][i];

                        found = true;
                        break;
                    }
                }

                // Create new pair by calculating xi and We

                if (!found)
                {
                    // New pair. Record xi and Weber.

                    pairsNew[c].p1 = i;
                    pairsNew[c].p2 = j;

                    #if dimension == 2
                    double ri = sqrt(V[i]/pi);
                    double rj = sqrt(V[j]/pi);
                    #else
                    double ri = cbrt(V[i]*3.0/pi/4.0);
                    double rj = cbrt(V[j]*3.0/pi/4.0);
                    #endif

                    pairsNew[c].xi = min(ri,rj)/max(ri,rj);

                    double req = 2*ri*rj/(ri+rj);
                    coord Ur = {U[i].x-U[j].x, U[i].x-U[j].x, U[i].x-U[j].x};
                    double UrSquare = sq(Ur.x)+sq(Ur.y)+sq(Ur.z);

                    // By definition, the first phase is the dispersed phase so
                    // we need rho1 here. If NOTENSION is set, then we require
                    // the SIGMA variable. Otherwise, we can use that of f.

                    #if NOTENSION
                    pairsNew[c].We = rho1*UrSquare*req/SIGMA;
                    #else
                    pairsNew[c].We = rho1*UrSquare*req/f.sigma;
                    #endif

                    pairsNew[c].A1 = A[i][j];
                    pairsNew[c].A2 = A[j][i];
                }

                c++;
            }
        }
    }

    // Move pointer and free memory

    free(pairs);
    pairs = pairsNew;
    Npairs = NpairsNew;

    // Clear up memory

    free(V);

    free(A[0]);
    free(A);
}


// Event called at the start of a simulation

event defaults (i = 0)
{
    // Remove f from interfaces and add a new duplicate of f as first interface.
    // The f field is kept as the sum of all interfaces. Also add a tag field.

    interfaces = NULL;
    tags = NULL;

    addInterface();
    addTag();

    scalar fnew = interfaces[0];

    foreach()
    {
        fnew[] = f[];
    }

    // Color numbers should not be averaged/interpolated in some way. They
    // should be injected, which preserves round numbers.

    #if TREE
    colors.restriction = restriction_tag;
    colors.prolongation = refine_injection;
    colors.refine = refine_injection;
    colors.dirty = true;
    #endif
}

// Event called before solving the vof equation

event vof (i++)
{
    // Create an array to store the number of particles in each void fraction
    // field

    int Nps[list_len(interfaces)];

    for (int i = 0; i < list_len(interfaces); i++)
    {
        Nps[i] = 0;
    }

    // Update tag fields

    int c = 0;

    NpOld = Np;
    Np = 0;

    int oldToNewNums[NpOld];

    for (int i = 0; i < NpOld; i++)
        oldToNewNums[i] = -1;

    for (scalar fi in interfaces)
    {
        // Remove very small object

        remove_droplets(f=fi);

        // Initialize g as a mask that is 1 for non-zero volume fractions

        scalar g[];

        foreach()
            g[] = fi[] > 0.0001;

        // Translate the mask field into tag numbers

        Nps[c] = tag(g);

        // Adjust the tag numbers so that they are larger than the tag numbers
        // belonging to previous fields

        if (c > 0)
            foreach()
                if (g[] > 0)
                    g[] += Np;

        // Transfer tag while keeping track of particle number changes

        transferTag(tags[c], g, oldToNewNums);

        // Update counts

        Np += Nps[c];
        c++;
    }

    checkOldToNewNums(oldToNewNums, NpOld, Np);

    // Reset particle interface number array

    setInterfaceNums(Nps);

    // Prevent coalescence

    preventCoalescence();

    // Update collision pairs

    updateCollisionPairs(oldToNewNums);

    // Merge pairs based in a physical criterion

    mergePairs();
}

// After advection, collect all interfaces into the f field. This is needed for
// the calculation of mixture density and viscosity

event properties (i++)
{
    collectInterfaces();
}