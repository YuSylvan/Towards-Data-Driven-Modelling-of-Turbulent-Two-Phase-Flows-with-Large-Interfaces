#include "navier-stokes/centered.h"
#include "snapshot.h"
#include "createArgs.h"
#define AVERAGING 1

//#define mu(f) (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "strainRateFix.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "timestep.h"

#define CASE 1

#if CASE == 1
//gas
#define RHO2 1
#define MU2 1.887E-05
#define SIGMA 0.07
//liquid
#define RHO1 10
#define MU1 1.887E-04

#elif CASE == 2

//gas
#define RHO2 1
#define MU2 1.887E-05
#define SIGMA 1
//liquid
#define RHO1 10
#define MU1 1.887E-04

#endif

#define GRAVITY 9.81
#define PressureGradient 0.5
 
#define RETAU 180.0
#define REBULK 880.0

#define HALFHEIGHT 0.05
//#define MU (1.0/RETAU)

// Domain dimensions need to be integer multiples of each other because there's
// no stretching.

#define HEIGHT (2.0*HALFHEIGHT)
#define LENGTH (4.0*HEIGHT)
#define WIDTH  (2.0*HEIGHT)

//#define HEIGHT (2.0*HALFHEIGHT)
//#define LENGTH (8.0*HEIGHT)
//#define WIDTH  (4.0*HEIGHT)

// Domain decomposition is proportional to the domain size, with cubic blocks
// for each processor. Decomposition is fully defined by the number of
// processors in the wall-normal direction.

#define NPROCY 2

// Global variables

face vector muv[];
face vector av[];
bool master;


int Nx, Ny, Nz;
double dx, dy, dz;
#if AVERAGING

    int Na;
    int Nr;
    int Nt;

    //scalar f[];

    // #define mu(f) MU
    // #define rho(f) RHO

    #define AVERAGEFREQ 10

    (const) scalar phis;

    // coord G = {0.0,0.0,0.0};
    // coord Z = {0.0,0.0,0.0};

    #include "twoPhaseAveraging.h"
    scalar fg[];
    scalar fl[];

    CREATEMEANFIELDS(fg)
    CREATEMEANFIELDS(fl)

#endif

int main (int argc, char * argv[])
{
    #include "setMaster.h"
    #include "readArgs.h"
    CFL = 0.3;
    a = av;
    L0 = LENGTH;
    periodic(right);

    #if dimension > 2
    periodic(front);
    #endif

    #if dimension == 2

    dimensions
    (
        nx=NPROCY*LENGTH/HEIGHT,
        ny=NPROCY
    );

    #else

    dimensions
    (
        nx=NPROCY*LENGTH/HEIGHT,
        ny=NPROCY,
        nz=NPROCY*WIDTH/HEIGHT
    );

    #endif

    origin(0.0, -HEIGHT/2.0, -WIDTH/2.0);
    N = pow(2,maxLevel);
    Nx = N;
    Ny = N*HEIGHT/LENGTH;
    Nz = N*WIDTH/LENGTH;

    dx = LENGTH/Nx;
    dy = HEIGHT/Ny;
    dz = WIDTH/Nz;

    rho1 = RHO1;
    rho2 = RHO2;

    mu1 = MU1;
    mu2 = MU2;

    f.sigma = SIGMA;

    //G.y = -GRAVITY;
    TOLERANCE = 1e-4;

    #if AVERAGING

    fg.nodump = true;
    fl.nodump = true;

    Na = N;
    Nr = N*HEIGHT/LENGTH;
    Nt = N*WIDTH/LENGTH;

    #endif

    run();
}



u.n[top] = dirichlet(0.0);
u.t[top] = dirichlet(0.0);

u.n[bottom] = dirichlet(0.0);
u.t[bottom] = dirichlet(0.0);



#if dimension > 2

u.r[top] = dirichlet(0.0);
u.r[bottom] = dirichlet(0.0);

#endif


#if AVERAGING

fg[left] = u.n[left] > 0.0 ? dirichlet(0.0) : neumann(0.0);
fg[right] = dirichlet(0.0);

fl[left] = u.n[left] > 0.0 ? dirichlet(1.0) : neumann(0.0);
fl[right] = dirichlet(1.0);

#endif




event properties (i++) 
{
  // fixme: metric
  foreach()
    rhov[] = rho(f[]);
  foreach_face () {
    alphav.x[] = 2./(rhov[] + rhov[-1]);
    double ff = (f[] + f[-1])/2.;
    muv.x[] = fm.x[]*mu(ff);
  }
}



event init (t = 0)
{

    if (!restore(file = "restart"))
    {
        foreach()
        {
            u.x[] =    1.5
              * (
                    3.0/2.0*(1.0-sq(y/HALFHEIGHT))
                  + sin(z*4.0*pi)/4.0
                );

            u.z[] = 1.5 * (sin(x*2.0*pi)/4.0);
        }

        fraction
        (
            f,
            -y

        );

        snapshot({f,p}, {u}, true, true, false);

    }

    #if AVERAGING

    foreach()
        // f[] = 1.0;

    resetAverages();

    #endif

    int nCells = 0;

    foreach(reduction(+:nCells))
    {
        nCells++;
    }

    fprintf(stdout, "Number of cells = %d\n", nCells);
}

#if AVERAGING

event average (i++)
{
    if (i%AVERAGEFREQ == 0)
    {
        double T = (t-tAvStart);

        if (T > 0)
        {
            double wa = (tAvPrev-tAvStart)/(T+1e-16);
            double wb = (t-tAvPrev+1e-16)/(T+1e-16);

            #include "twoPhaseAveragingUpdate.h"
            foreach()
            {
                fg[] = clamp(1.0-f[],0,1);
                fl[] = clamp(f[],0,1);
            }

            PHASEAVERAGE(fl)
            PHASEAVERAGE(fg)

            tAvPrev = t;
        }
    }
}

#endif

event velocity (i++,last) {
  dt = dtnext (timestep (u, 1));
}
// event stability (i++)
// {
//     CFL = 0.3;

// }

event acceleration (i++, last) 

{
    
    foreach_face(x)
    {
        av.x[] += PressureGradient / rho[];
    }

    foreach_face(y)
    {  
       av.y[] -= GRAVITY;
    }
        
}


event logfile (i++)
{

    double A = 0.0;
    double Ub = 0.0;

    foreach_face
    (
        x,
        reduction(+:A)
        reduction(+:Ub)
    )
    {
        if
        (
            x > L0/2.0-Delta/2.0
         && x < L0/2.0+Delta/2.0
        )
        {
            double Ai = sq(Delta)*fm.x[];

            A += Ai;
            Ub += face_value(u.x,0)*Ai;
        }
    }

    Ub /= A;

    if (master)
    {
        fprintf
        (
            stdout,
            "i = %d, t = %.8g, dt = %.8g, Ub = %.8g\n",
            i,
            t,
            dt,
            Ub
        );

        fflush(stdout);
    }
}



event snapshots (t = writeTime; t <= endTime; t += writeTime)
{
    scalar P[];
    // Save snapshot with primitive variables
    foreach()
    {
        P[] = p[] + rho(f[])*(G.x*(x-Z.x) + G.y*(y-Z.y) + G.z*(z-Z.z));
    }

    #if AVERAGING

        SETNODUMPMEANFIELDS(fg)
        SETNODUMPMEANFIELDS(fl)

        const int dirIndex = snapshot({f,p,P,rho}, {u}, true, true, true);

        #include "twoPhaseAveragingCartesian.h"

        PHASESNAPSHOT_CARTESIAN(fl)
        PHASESNAPSHOT_CARTESIAN(fg)

        resetAverages();

    #else

        snapshot({f,p,rho}, {u}, true, true, true);

    #endif
}

#define outputtime 0.02

event output (t = outputtime; t <= endTime; t += outputtime)

{
  scalar omega[];
  vorticity (u, omega);
  char filename[25] = {"./Line/X"};
  char timeindex[10];
  sprintf(timeindex,"%.3lf",t);
  strcat(filename, timeindex);
  strcat(filename, "s.txt");
  FILE * fp = fopen (filename, "w");
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({u,f,p,omega }, fp, n = N, linear = true,box = {{0.,-0.05},{L0,0.05}});
}

event output2 (t = outputtime; t <= endTime; t += outputtime)
{
  char filename[25] = {"./Line/Y"};
  char timeindex[10];
  sprintf(timeindex,"%.3lf",t);
  strcat(filename, timeindex);
  strcat(filename, "s.txt");
  FILE * fp = fopen (filename, "w");
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({u,f,p}, fp, n = 5, linear = true,box = {{0.5001,-0.05},{0.5060,0.05}});
}