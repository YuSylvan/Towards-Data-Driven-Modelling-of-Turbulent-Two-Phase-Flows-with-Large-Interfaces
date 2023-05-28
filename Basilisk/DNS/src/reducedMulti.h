// The default reduced.h only works on f. However, when having multiple
// interfaces for the same phase, f is not part of the interfaces array and
// buoyancy is thus not taken into account. We need to add the gravity potential
// to each interface separately.

coord G = {0.,0.,0.}, Z = {0.,0.,0.};

#include "iforce.h"
#include "curvature.h"

event acceleration (i++)
{
    coord G1;

    foreach_dimension()
        G1.x = (rho2 - rho1)*G.x;

    for (scalar f in interfaces)
    {
        scalar phi = f.phi;

        if (phi.i)
        {
            position (f, phi, G1, Z, add = true);
        }
        else
        {
            phi = new scalar;
            position (f, phi, G1, Z, add = false);
            f.phi = phi;
        }
    }
}
