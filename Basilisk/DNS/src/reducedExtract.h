// The buoyancy potential is added to f.phi. In order to extract it, we need to
// compute the difference in phi between its values before and after the
// buoyancy acceleration event. We assume that there is only one interface (f).
// Note that execution of events seems to be in reversed order (check with "qcc
// -events ...").

scalar phig[];

// By including iforce.h here, we assure that its acceleration event is executed
// last. This is important.

#include "iforce.h"

event acceleration (i++, last)
{
    // This event is executed last. Extract the phi value.

    scalar phi = f.phi;

    foreach()
        phig[] += phi[];
}

#include "reduced.h"

event acceleration (i++)
{
    // This event is executed first. Store the current value of phi.

    scalar phi = f.phi;

    if (phi.i)
    {
        foreach()
            phig[] = -phi[];
    }
    else
    {
        foreach()
            phig[] = 0.0;
    }
}
