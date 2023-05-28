// This solver has three modes of operation, set by the epsilon variable (third
// argument):

//  (1) Start from a uniform velocity field and initializing a bubble
//  (2) Resume from a single phase turbulent field generated by the pipeFlow app
//  (3) Resume run

#define mu(f) (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#define REDUCEDRUN 0

#include "embed.h"
#include "navier-stokes/centered.h"
#include "strainRateFix.h"
#include "two-phase.h"
#include "reduced.h"
#include "tensionExtract.h"
#include "snapshot.h"
#include "gradingFactor.h"
#include "createArgs.h"
#include "tag.h"

// Non-dimensional numbers are defined as RE = RHO1*U*D/MU1, WE =
// RHO1*U^2*D/SIGMA, FR = U/sqrt(GRAVITY*D), GAMMARHO = RHO1/RHO2 and GAMMAMU =
// MU1/MU2. These are five non-dimensional numbers (RE, WE, FR, GAMMARHO,
// GAMMAMU) with eight dimensional parameters (RHO1, RHO2, MU1, MU2, U, D,
// SIGMA, GRAVITY), so we can set three dimensional parameters to unity. We
// conveniently hard-code RHO1 = 1, U = 1 and D = 1, so that the non-dimensional
// numbers reduce to: RE = 1/MU1, WE = 1/SIGMA, FR = 1/sqrt(GRAVITY), GAMMAMU =
// MU1/MU2 and GAMMARHO = 1/RHO2.

#define RE 5300.0
#define WE 40.0
#define FR 1.5
#define GAMMARHO 10.0
#define GAMMAMU 10.0

// Five non-dimensional numbers give five dimensional parameters

#define MU1 (1.0/RE)
#define MU2 (MU1/GAMMAMU)
#define SIGMA (1.0/WE)
#define GRAVITY (1.0/(FR*FR))
#define RHO2 (1.0/GAMMARHO)

// The co- or counter-current liquid velocity can be set here in units of U. For
// UL = 1 or -1, the single phase liquid flow has Reynolds number equal to RE.

#define UL 1.0

#define RECYCLING 1
#define LENGTH 16.0

#define SAMPLEDIST 4.0

#define LB0 2.0
#define DB0 0.85
#define XB0 (LENGTH-SAMPLEDIST-LB0/2.0-2.0)

bool master;

// Bubble velocity w.r.t. the inertial frame, inlet velocity and wall velocity

double Ub = 0.0;
double Ui = 0.0;
double Uw = 0.0;

// Helper field for the inlet BC

vector v[];

#if !REDUCEDRUN

    int Na;
    int Nr;
    int Nt;

    #define AVERAGEFREQ 10
    #define BUBBLESTATSFREQ 100

    #include "twoPhaseAveraging.h"

    scalar fg[];
    scalar fl[];

    CREATEMEANFIELDS(fg)
    CREATEMEANFIELDS(fl)

#endif

#include "functions.c"

int main (int argc, char * argv[])
{
    #include "setMaster.h"
    #include "readArgs.h"

    foreach_dimension()
    {
        v.x.nodump = true;
    }

    #if !REDUCEDRUN
    fg.nodump = true;
    fl.nodump = true;
    #endif

    dimensions(nx=8, ny=1, nz=1);

    size(LENGTH);
    origin(0.0, -0.5, -0.5);

    init_grid(pow(2,maxLevel));

    // Physical parameters

    rho1 = 1.0;
    rho2 = RHO2;

    mu1 = MU1;
    mu2 = MU2;

    f.sigma = SIGMA;

    G.x = -GRAVITY;
    Z.x = XB0;

    // Cylindrical grid parameters

    #if !REDUCEDRUN
    Na = pow(2, maxLevel);
    Nr = ceil(0.5/L0*Na)*2;
    Nt = ceil(pi*Nr);
    #endif

    NITERMAX = 20;

    run();
}

// For recycling, use the helper field v. The helper field v can be updated
// using the recycle() function. If no recycling is used, a plug flow boundary
// condition is directly imposed.

f[left] = u.n[left] > 0.0 ? dirichlet(1.0) : neumann(0.0);
f[right] = dirichlet(1.0);

#if !REDUCEDRUN

fg[left] = u.n[left] > 0.0 ? dirichlet(0.0) : neumann(0.0);
fg[right] = dirichlet(0.0);

fl[left] = u.n[left] > 0.0 ? dirichlet(1.0) : neumann(0.0);
fl[right] = dirichlet(1.0);

#endif

p[right] = neumann(0.0);
p[left] = dirichlet(0.0);

// On the left boundary, limit the inflow to zero such that the left boundary
// always remains an outflow boundary. This is necessary to prevent
// instabilities. It is only physical if the liquid bulk velocity is negative.
// In single phase flow, this is guaranteed by imposing a fictitious bubble
// velocity of 0.4 (see the bubbleVelocity() function). In multiphase flow, the
// bubble rise velocity is always larger than the bulk velocity, so that in the
// MFR the bulk velocity is indeed negative.

u.n[right] = dirichlet(RECYCLING ? v.x[] : Ui);
u.n[embed] = dirichlet(Uw);
u.n[left] = u.n[] > 0.0 ? dirichlet(0.0) : neumann(0.0);
u.n[top] = dirichlet(0.0);
u.n[bottom] = dirichlet(0.0);

u.t[right] = dirichlet(RECYCLING ? v.y[] : 0.0);
u.t[embed] = dirichlet(0.0);
u.t[left] = neumann(0.0);

// Table that gives the relation between (n,t) and (x,y) in 2D:
//
//    | left  top
//    | right bottom
// -----------------------
//  n | x     y
//  t | y     x

// Table that gives the relation between (n,t,r) and (x,y,z) n 3D:
//
//    | left  top    front
//    | right bottom back
// ------------------------------
//  n | x     y      z
//  t | y     z      x
//  r | z     x      y

u.t[top] = dirichlet(0.0);
u.t[bottom] = dirichlet(0.0);

u.n[front] = dirichlet(0.0);
u.n[back] = dirichlet(0.0);

u.t[front] = dirichlet(Uw);
u.t[back] = dirichlet(Uw);

u.r[right] = dirichlet(RECYCLING ? v.z[] : 0.0);
u.r[embed] = dirichlet(0.0);
u.r[left] = neumann(0.0);
u.r[top] = dirichlet(Uw);
u.r[bottom] = dirichlet(Uw);
u.r[front] = dirichlet(0.0);
u.r[back] = dirichlet(0.0);

event init(t = 0)
{
    if (epsilon == 1)
    {
        // Start from a uniform velocity field and initializing a bubble

        if (master)
        {
            fprintf(stdout, "Multiphase run, starting from plug flow\n");
        }

        updateApertures();
        initLiquid();
        updateFlux();
        initPlugVelocity(0.0);
        updateFlux();
        initBubble();
        updateFlux();

        // The bubble is initialized with zero velocity, so w.r.t. the inertial
        // frame it is moving at -Uw

        Ui = 0.0;
        Uw = -UL;
        Ub = -Uw;

        snapshot({cs, f}, {u}, true, true, false);
    }
    else if (epsilon == 2)
    {
        // Resume from a turbulent single phase field generated by the pipeFlow
        // app

        if (master)
        {
            fprintf(stdout, "Multiphase run, resuming from a turbulent single phase field\n");
        }

        if (!restore(file = "restart"))
        {
            if (master)
            {
                fprintf(stdout, "Could not find restart field\n");
            }

            exit(1);
        }

        // The pipeFlow app simulates on a different domain. Rescale and shift
        // the domain.

        size(LENGTH);
        origin(0.0, -0.5, -0.5);

        updateApertures();
        initLiquid();
        updateFlux();

        // Scale the velocity by its bulk velocity, bringing the velocity scale
        // back to unity. Then shift so that its streamwise mean is zero.

        scaleAndShiftVelocity(liquidVelocity(L0/2.0));
        updateFlux();

        // Initialize the bubble. This will destroy continuity of the turbulent
        // field and thus induce an initial error, which should gradually decay.

        initBubble();
        updateFlux();

        // The bubble is initialized with zero velocity, so w.r.t. the inertial
        // frame it is moving at -Uw

        Ui = 0.0;
        Uw = -UL;
        Ub = -Uw;
    }
    else if (epsilon == 3)
    {
        // Resume run

        if (master)
        {
            fprintf(stdout, "Resume run\n");
        }

        if (!restore(file = "restart"))
        {
            if (master)
            {
                fprintf(stdout, "Could not find restart field\n");
            }

            exit(1);
        }

        updateApertures();
        updateFlux();

        Ui = liquidVelocity(X0+L0-SAMPLEDIST);
        Uw = Ui-UL;
        Ub = bubbleVelocity(XB0+LB0/2.0)-Uw;
    }
    else
    {
        if (master)
        {
            fprintf(stdout, "Invalid value for epsilon (should be 1, 2 or 3)");
        }

        exit(1);
    }

    #if RECYCLING
    generateRecyclePoints();
    #endif

    int nCells = 0;

    foreach(reduction(+:nCells))
    {
        nCells++;
    }

    fprintf(stdout, "Number of cells = %d\n", nCells);

    #if !REDUCEDRUN
    resetAverages();
    #endif
}

event stability (i++)
{
    CFL = 0.4;
}

// At the beginning of the time step, adjust the MFR

double Uslip = 0;

event mfr (i++)
{
    Uslip = bubbleVelocity(XB0+LB0/2.0);

    Ui -= Uslip;
    Uw = Ui-UL;
    Ub = -Uw;

    #if RECYCLING
    recycle();
    #endif

    boundary((scalar*){u});
}

// MFR acceleration implementation of the temporal term

event acceleration (i++)
{
    foreach_face(x)
        a.x[] -= Uslip/dt;
}

// MFR acceleration implementation of the spatial term, implemented by
// translation of the face flux into the inertial one

event advection_term (i++)
{
    // Make flux relative to the inertial frame

    foreach_face(x)
        uf.x[] += fm.x[]*Ub;
}

event viscous_term (i++)
{
    // Make flux relative to the moving frame

    foreach_face(x)
        uf.x[] -= fm.x[]*Ub;
}

#if !REDUCEDRUN

event average (i++)
{
    // Averaging

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

    // Bubble statistics

    if (i%BUBBLESTATSFREQ == 0)
    {
        coord Pb;
        double Vb = bubbleVolume(&Pb);

        if (master)
        {
            fprintf
            (
                stdout,
                "i = %6d, t = %12.8f, Vb = %14.8e, xb = %14.8e, yb = %15.8e, zb = %15.8e\n",
                i, t, Vb, Pb.x, Pb.y, Pb.z
            );
        }
    }
}

#endif

// Runtime output

int nCells = -1;

event logfile (i++)
{
    double Ul = liquidVelocity(X0+L0-SAMPLEDIST);

    if (master)
    {
        fprintf
        (
            stdout,
            "i = %6d, t = %12.8f, dt = %14.8e, Ub = %14.8e, Ui = %14.8e, Ul = %14.8e, Uslip = %14.8e\n",
            i,
            t,
            dt,
            Ub,
            Ui,
            Ul,
            Uslip
        );

        fflush(stdout);
    }
}

// Take snapshot

event snapshots (t = writeTime; t <= endTime; t += writeTime)
{
    // Calculate absolute pressure and store

    scalar P[];

    foreach()
    {
        P[] = p[] + rho(f[])*(G.x*(x-Z.x) + G.y*(y-Z.y) + G.z*(z-Z.z));
    }

    // Save snapshot with primitive variables

    #if !REDUCEDRUN

        SETNODUMPMEANFIELDS(fg)
        SETNODUMPMEANFIELDS(fl)

        const int dirIndex = snapshot({f,P}, {u}, true, true, true);

        #include "twoPhaseAveragingCylindrical.h"

        PHASESNAPSHOT_CYLINDRICAL(fl)
        PHASESNAPSHOT_CYLINDRICAL(fg)

        resetAverages();

    #else

        snapshot({f,P}, {u}, true, true, true);

    #endif
}