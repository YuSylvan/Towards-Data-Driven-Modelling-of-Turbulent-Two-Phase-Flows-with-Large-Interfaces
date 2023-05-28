// Compute required fields

scalar r[];
vector surf[];

foreach()
{
    r[] = sqrt(sq(y)+sq(z));
}

if (!is_constant(phis))
{
    face vector as[];

    reconstructAcceleration(as,phis,f);
    reconstruct(surf,as);
}
else
{
    foreach()
    {
        foreach_dimension()
        {
            surf.x[] = 0.0;
        }
    }
}

face vector nu[];

foreach_face()
{
    nu.x[] = mu.x[]/rho(face_value(f,0));
}

vector dnuduidjdj[];

lapl(dnuduidjdj.x, nu, u.x);
lapl(dnuduidjdj.y, nu, u.y);
lapl(dnuduidjdj.z, nu, u.z);

vector duxdj[];
vector duydj[];
vector duzdj[];

grad(duxdj,u.x);
grad(duydj,u.y);
grad(duzdj,u.z);