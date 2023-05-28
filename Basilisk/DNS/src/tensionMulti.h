// The default tension.h stability event limits the time step by the sum of all
// surface tension values of all interfaces. When using multiple interfaces for
// the same phase, it is sufficient to limit the time step by the actual surface
// tension value sigma. We add an event that reduces the surface tension such
// that the average across all interfaces is sigma, and then resets the value.
// Note that execution of events seems to be in reversed order (check with "qcc
// -events ...").

event stability (i++)
{
    // This event is executed last. Scale back the sigma value.

    for (scalar c in interfaces)
        c.sigma = c.sigma*list_len(interfaces);
}

#include "tension.h"

event stability (i++)
{
    // This event is executed first. Scale the sigma value by the number of
    // interfaces.

    for (scalar c in interfaces)
        c.sigma = c.sigma/list_len(interfaces);
}