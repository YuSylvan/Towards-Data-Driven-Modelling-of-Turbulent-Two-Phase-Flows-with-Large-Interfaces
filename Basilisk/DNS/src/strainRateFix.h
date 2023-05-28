// In case of embedding, the strain rate tensor is approximated as Sij =
// 1/2*dui/dj. When having variable viscosity, a term 1/rho*dmu/dj*duj/di is
// then missing from the momentum equation. This term can be added here
// explicitly, which is only useful in case of variable viscosity.

// This file should be included after navier-stokes/centered.h, so that the
// event is executed at the right time.

#if EMBED

event viscous_term (i++)
{
    if (!is_constant(mu.x))
    {
        correction(dt);

        foreach()
        {
            foreach_dimension()
            {

                u.x[] +=
                    dt/(2.0*cm[]*rho[]*sq(Delta)+SEPS)
                  * (
                        (mu.x[1] - mu.x[])*(fm.x[1]*u.x[1] - fm.x[]*u.x[-1])
                        #if dimension > 1
                      + (mu.y[0,1] - mu.y[])*(fm.x[1]*u.y[1] - fm.x[]*u.y[-1])
                        #endif
                        #if dimension > 2
                      + (mu.z[0,0,1] - mu.z[])*(fm.x[1]*u.z[1] - fm.x[]*u.z[-1])
                        #endif
                    );
            }
        }

        correction(-dt);
    }
}

#endif