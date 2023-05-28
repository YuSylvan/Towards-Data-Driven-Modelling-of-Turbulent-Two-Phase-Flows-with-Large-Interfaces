// Based on http://basilisk.fr/sandbox/weugene/output_vtu_foreach.h

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

// For some obscure reason, this file needs to be included like this

@include <dirent.h>

#ifndef DEFORMMESH

#define DEFORMX(x) (x)
#define DEFORMY(y) (y)
#define DEFORMZ(z) (z)

#endif

// Write the vertex index and increase the counter

#define WRITE_INDEX(i,j,k) fprintf(fp, "%d ", l++);

// Compute coordinate adjusted for periodicity

#if MULTIGRID_MPI

    #if dimension == 2

    #define LX0 (L0*(double)(mpi_dims[0])/max(mpi_dims[0],mpi_dims[1]))
    #define LY0 (L0*(double)(mpi_dims[1])/max(mpi_dims[0],mpi_dims[1]))
    #define LZ0 L0

    #else

    #define LX0 (L0*(double)(mpi_dims[0])/max(max(mpi_dims[0],mpi_dims[1]),mpi_dims[2]))
    #define LY0 (L0*(double)(mpi_dims[1])/max(max(mpi_dims[0],mpi_dims[1]),mpi_dims[2]))
    #define LZ0 (L0*(double)(mpi_dims[2])/max(max(mpi_dims[0],mpi_dims[1]),mpi_dims[2]))

    #endif

#else

    #define LX0 L0;
    #define LY0 L0;
    #define LZ0 L0;

#endif

#define PERIODICX(i,j,k,xijk) (i == 1 && xijk < x ? DEFORMX(xijk+LX0) : DEFORMX(xijk))
#define PERIODICY(i,j,k,yijk) (j == 1 && yijk < y ? DEFORMY(yijk+LY0) : DEFORMY(yijk))
#define PERIODICZ(i,j,k,zijk) (k == 1 && zijk < z ? DEFORMZ(zijk+LZ0) : DEFORMZ(zijk))

// Write vertex coordinates

#define WRITE_VERTEX(i,j,k)                     \
{                                               \
    double xijk = xv[i,j,k];                    \
    double yijk = yv[i,j,k];                    \
    double zijk = zv[i,j,k];                    \
                                                \
    fprintf(fp, "%g ", PERIODICX(i,j,k,xijk));  \
    fprintf(fp, "%g ", PERIODICY(i,j,k,yijk));  \
    fprintf(fp, "%g ", PERIODICZ(i,j,k,zijk));  \
    fprintf(fp, "\n");                          \
}

// Write binary vertex coordinates

#define WRITE_VERTEX_BINARY(i,j,k)              \
{                                               \
    double xijk = xv[i,j,k];                    \
    double yijk = yv[i,j,k];                    \
    double zijk = zv[i,j,k];                    \
                                                \
    double xi = PERIODICX(i,j,k,xijk);          \
    double yj = PERIODICY(i,j,k,yijk);          \
    double zk = PERIODICZ(i,j,k,zijk);          \
                                                \
    fwrite(&xi, sizeof(double), 1, fp);         \
    fwrite(&yj, sizeof(double), 1, fp);         \
    fwrite(&zk, sizeof(double), 1, fp);         \
}

void output_pvtu
(
    scalar * list,
    vector * vlist,
    FILE * fp,
    bool binary
)
{
    fputs("<?xml version=\"1.0\"?>\n", fp);
    fputs("<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

    fputs("  <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs("      <PCellData Scalars=\"scalars\">\n", fp);

    char format[80];

    if (binary)
    {
        sprintf(format, "append");
    }
    else
    {
        sprintf(format, "ascii");
    }

    for (scalar s in list)
    {
        fprintf(fp, "        <PDataArray type=\"Float64\" Name=\"%s\" format=\"%s\">\n", s.name, format);
        fputs("        </PDataArray>\n", fp);
    }

    for (vector v in vlist)
    {
        char vName[strlen(v.x.name)];
        strcpy(vName, v.x.name);
        vName[strlen(v.x.name)-2] = '\0';

        fprintf(fp, "        <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"%s\">\n", vName, format);
        fputs("        </PDataArray>\n", fp);
    }

    fputs("      </PCellData>\n", fp);

    fputs("      <PPoints>\n", fp);
    fputs("        <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs("        </PDataArray>\n", fp);
    fputs("      </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
    {
        fprintf(fp, "      <Piece Source=\"data_%04d.vtu\"/>\n", i);
    }

    fputs("  </PUnstructuredGrid>\n", fp);

    fputs("</VTKFile>\n", fp);
}

void output_vtu_foreach
(
    scalar * list,
    vector * vlist,
    int n,
    FILE * fp,
    bool binary
)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif

    vertex scalar xv[];
    vertex scalar yv[];
    vertex scalar zv[];

    foreach_vertex()
    {
        xv[] = x;
        yv[] = y;
        zv[] = z;
    }

    int nCells = 0;

    foreach(serial)
    {
        nCells++;
    }

    #if dimension == 1
    const int nPoints = 2*nCells;
    #elif dimension == 2
    const int nPoints = 4*nCells;
    #else
    const int nPoints = 8*nCells;
    #endif

    fputs("<?xml version=\"1.0\"?>\n", fp);
    fputs("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs("  <UnstructuredGrid>\n", fp);

    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nPoints, nCells);
    fputs("      <CellData Scalars=\"scalars\">\n", fp);

    // Scalar data

    int byteCount = 0;

    for (scalar s in list)
    {
        if (binary)
        {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name, byteCount);
            fputs("        </DataArray>\n", fp);

            byteCount += (nCells+1)*8;
        }
        else
        {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);

            foreach()
            {
                fprintf(fp, "%g \n", val(s));
            }

            fputs("        </DataArray>\n", fp);
        }
    }

    // Vector data

    for (vector v in vlist)
    {
        char vName[strlen(v.x.name)];
        strcpy(vName, v.x.name);
        vName[strlen(v.x.name)-2] = '\0';

        if (binary)
        {
            fprintf(fp,"        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", vName, byteCount);
            fputs("        </DataArray>\n", fp);

            byteCount += (nCells*3+1)*8;
        }
        else
        {
            fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", vName);

            foreach()
            {
                #if dimension == 1
                fprintf(fp, "%g 0.0 0.0 \n", val(v.x));
                #elif dimension == 2
                fprintf(fp, "%g %g 0.0 \n", val(v.x), val(v.y));
                #elif dimension == 3
                fprintf(fp, "%g %g %g \n", val(v.x), val(v.y), val(v.z));
                #endif
            }

            fputs("        </DataArray>\n", fp);
        }
    }

    fputs("      </CellData>\n", fp);

    // Point data

    fputs("      <Points>\n", fp);

    if (binary)
    {
        fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", byteCount);

        byteCount += (nPoints*3+1)*8;

        fputs("        </DataArray>\n", fp);
    }
    else
    {
        fputs("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

        foreach()
        {
            WRITE_VERTEX(0,0,0)
            WRITE_VERTEX(1,0,0)

            #if dimension > 1
            WRITE_VERTEX(1,1,0)
            WRITE_VERTEX(0,1,0)
            #endif

            #if dimension > 2
            WRITE_VERTEX(0,0,1)
            WRITE_VERTEX(1,0,1)
            WRITE_VERTEX(1,1,1)
            WRITE_VERTEX(0,1,1)
            #endif
        }

        fputs("        </DataArray>\n", fp);
    }

    fputs("      </Points>\n", fp);

    // Cell connectivity

    fputs("      <Cells>\n", fp);
    fputs("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);

    int l = 0;

    foreach(serial)
    {
        WRITE_INDEX(0,0,0)
        WRITE_INDEX(1,0,0)

        #if dimension > 1
        WRITE_INDEX(1,1,0)
        WRITE_INDEX(0,1,0)
        #endif

        #if dimension > 2
        WRITE_INDEX(0,0,1)
        WRITE_INDEX(1,0,1)
        WRITE_INDEX(1,1,1)
        WRITE_INDEX(0,1,1)
        #endif

        fprintf(fp, "\n");
    }

    fputs("        </DataArray>\n", fp);

    // Offsets

    fputs("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

    for (int i = 1; i < nCells+1; i++)
    {
        #if dimension == 1
        fprintf(fp, "%d \n", i*2);
        #elif dimension == 2
        fprintf(fp, "%d \n", i*4);
        #else
        fprintf(fp, "%d \n", i*8);
        #endif
    }

    fputs("        </DataArray>\n", fp);

    // Cell types

    fputs("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);

    foreach()
    {
        #if dimension == 1
        fputs("3\n", fp);
        #elif dimension == 2
        fputs("9\n", fp);
        #else
        fputs("12\n", fp);
        #endif
    }

    fputs("        </DataArray>\n", fp);
    fputs("      </Cells>\n", fp);
    fputs("    </Piece>\n", fp);
    fputs("  </UnstructuredGrid>\n", fp);

    if (binary)
    {
        fputs("  <AppendedData encoding=\"raw\">\n", fp);
        fputs("_", fp);

        unsigned long long blockLength = nCells*8;

        #if dimension < 2
        double vy = 0.0;
        #endif

        #if dimension < 3
        double vz = 0.0;
        #endif

        // Binary scalar data

        for (scalar s in list)
        {
            fwrite(&blockLength, sizeof(unsigned long long), 1, fp);

            foreach()
            {
                fwrite(&val(s), sizeof(double), 1, fp);
            }
        }

        // Binary vector data

        blockLength = nCells*8*3;

        for (vector v in vlist)
        {
            fwrite(&blockLength, sizeof(unsigned long long), 1, fp);

            foreach()
            {
                #if dimension == 1
                fwrite(&val(v.x), sizeof(double), 1, fp);
                fwrite(&vy, sizeof(double), 1, fp);
                fwrite(&vz, sizeof(double), 1, fp);
                #elif dimension == 2
                fwrite(&val(v.x), sizeof(double), 1, fp);
                fwrite(&val(v.y), sizeof(double), 1, fp);
                fwrite(&vz, sizeof(double), 1, fp);
                #else
                fwrite(&val(v.x), sizeof(double), 1, fp);
                fwrite(&val(v.y), sizeof(double), 1, fp);
                fwrite(&val(v.z), sizeof(double), 1, fp);
                #endif
            }
        }

        // Binary point data

        blockLength=nPoints*8*3;

        fwrite(&blockLength, sizeof(unsigned long long), 1, fp);

        foreach()
        {
            WRITE_VERTEX_BINARY(0,0,0)
            WRITE_VERTEX_BINARY(1,0,0)

            #if dimension > 1
            WRITE_VERTEX_BINARY(1,1,0)
            WRITE_VERTEX_BINARY(0,1,0)
            #endif

            #if dimension > 2
            WRITE_VERTEX_BINARY(0,0,1)
            WRITE_VERTEX_BINARY(1,0,1)
            WRITE_VERTEX_BINARY(1,1,1)
            WRITE_VERTEX_BINARY(0,1,1)
            #endif
        }

        fputs("\n", fp);
        fputs("  </AppendedData>\n", fp);
    }

    fputs("</VTKFile>\n", fp);

    fflush(fp);

    #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
    #endif
}

const char PvdFileName[] = "data.pvd";

void readPvd(Array* indices, Array* times)
{
    if (access(PvdFileName, F_OK) != -1)
    {
        FILE* PvdFile;

        PvdFile = fopen(PvdFileName, "r");

        int i;
        double time;

        char buffer[1024];
        char * ptr;

        while ((ptr = fgets(buffer, sizeof(buffer), PvdFile)) != NULL)
        {
            if (sscanf(buffer, "    <DataSet timestep=\"%lf\" file=\"%d/data.pvtu\"/>\n", &time, &i) == 2)
            {
                array_append(times, &time, sizeof(double));
                array_append(indices, &i, sizeof(int));
            }
        }

        fclose(PvdFile);
    }
}

void writePvd(Array* indices, Array* times)
{
    FILE * PvdFile;

    PvdFile = fopen(PvdFileName, "w");

    fputs("<?xml version=\"1.0\"?>\n", PvdFile);
    fputs("<VTKFile type=\"Collection\" version=\"0.1\">\n", PvdFile);
    fputs("  <Collection>\n", PvdFile);

    if (times->len)
    {
        const double * ts = (double*) times->p;
        const int * is = (int*) indices->p;

        const int numLines = times->len/sizeof(double);

        for (int j = 0; j < numLines; j++)
        {
            const int i = is[j];
            const double time = ts[j];

            fprintf(PvdFile, "    <DataSet timestep=\"%g\" file=\"%06d/data.pvtu\"/>\n", time, i);
        }
    }

    fputs("  </Collection>\n", PvdFile);
    fputs("</VTKFile>\n", PvdFile);

    fclose(PvdFile);
}

int snapshot
(
    scalar * slist,
    vector * vlist,
    bool binary,
    bool storeVTU,
    bool storeDump
)
{
    bool master;

    #if !_MPI
    master = true;
    #else
    master = (mpi_rank == 0);
    #endif

    // Read stored times

    Array* indices = array_new();
    Array* times = array_new();

    if (master)
    {
        readPvd(indices, times);
    }

    // Check if this time is already stored

    int alreadyStored = 0;

    if (times->len && master)
    {
        double * ts = (double*) times->p;

        const int numLines = times->len/sizeof(double);

        for (int i = 0; i < numLines; i++)
        {
            if (round(ts[i]/1e-8) == round(t/1e-8))
            {
                alreadyStored = true;
                break;
            }
        }
    }

    #if _MPI
    mpi_all_reduce(alreadyStored, MPI_INT, MPI_MAX);
    #endif

    int dirNum = -1;

    if (!alreadyStored)
    {
        if (master && times->len)
        {
            dirNum = times->len/sizeof(double);
        }
        else
        {
            dirNum = 0;
        }

        #if _MPI
        MPI_Bcast(&dirNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
        #endif

        char dirName[20];
        sprintf(dirName, "%06d", dirNum);

        if (master)
        {
            // Make directory

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
        }
        else
        {
            // Wait until the directory appears or sleep for 100 ms. This seems
            // to be necessary for consistency on network file systems

            DIR * dir;
            dir = NULL;

            struct timespec t;

            t.tv_sec = 0;
            t.tv_nsec = 0.1*1e9;

            while ((dir = opendir(dirName)) == NULL)
            {
                nanosleep(&t, NULL);
            }

            if (dir)
            {
                closedir(dir);
            }
        }

        #if _MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        if (storeDump)
        {
            // Store dump data

            char fileName[80];

            sprintf(fileName, "%s/dump", dirName);

            /*
                // Debug code to check which fields are being dumped:

                scalar * list = dump_list (all);
                for (scalar s in list)
                {
                    fprintf(stdout, "%s\n", s.name);
                }
            */

            dump(file = fileName);
        }

        if (storeVTU)
        {
            // Store VTU data

            FILE* fileVTU;

            char fileName[80];

            sprintf(fileName, "%s/data_%04d.vtu", dirName, pid());

            fileVTU = fopen(fileName, "w");

            output_vtu_foreach(slist, vlist, N, fileVTU, binary);

            fclose(fileVTU);

            if (master)
            {
                sprintf(fileName, "%s/data.pvtu", dirName);

                FILE* filePVTU;
                filePVTU = fopen(fileName, "w");

                output_pvtu(slist, vlist, filePVTU, binary);

                fclose(filePVTU);
            }
        }

        if (master)
        {
            // Update PVD file

            array_append(times, &t, sizeof(double));
            array_append(indices, &dirNum, sizeof(int));

            writePvd(indices, times);
        }

        #if _MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        if (master)
        {
            fprintf(stdout, "Saved snapshot for time %.8g to %s\n", t, dirName);
            fflush(stdout);
        }
    }
    else
    {
        if (master)
        {
            fprintf(stdout, "Tried to save snapshot for time %.8g but it is already stored\n", t);
        }
    }

    array_free(indices);
    array_free(times);

    return dirNum;
}

void writeSpatialDataAscii
(
    double ** coords,
    double * fields,
    const char fileName[],
    const int dirIndex,
    const int Nc,
    const int Nf,
    const int Ni,
    const bool writeCoords
)
{
    bool master;

    #if !_MPI
    master = true;
    #else
    master = (mpi_rank == 0);
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

    for (int i = 0; i < Ni; i++)
    {
        if (writeCoords)
        {
            for (int c = 0; c < Nc; c++)
            {
                fprintf(file, "%g ", coords[c][i]);
            }
        }

        for (int l = 0; l < Nf; l++)
        {
            fprintf(file, "%g ", fields[l*Ni+i]);
        }

        fprintf(file, "\n");
    }

    fclose(file);

    if (master)
    {
        fprintf(stdout, "Wrote data to %s\n", filePath);
        fflush(stdout);
    }
}
