// Macro to compute a vector transform

#define TRANSVEC(A,R,T,X,Y,Z)                                                   \
                                                                                \
foreach()                                                                       \
{                                                                               \
    double sin = y/sqrt(sq(y) + sq(z));                                         \
    double cos = z/sqrt(sq(y) + sq(z));                                         \
                                                                                \
    A[] = X[];                                                                  \
    R[] = sin*Y[] + cos*Z[];                                                    \
    T[] = cos*Y[] - sin*Z[];                                                    \
}

// Macro to compute a symmetric tensor transform

#define TRANSTENS(AA,AR,AT,RA,RR,RT,TA,TR,TT,XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ)        \
                                                                                \
foreach()                                                                       \
{                                                                               \
    double sin = y/sqrt(sq(y) + sq(z));                                         \
    double cos = z/sqrt(sq(y) + sq(z));                                         \
                                                                                \
    AA[] = XX[];                                                                \
    AR[] = sin*XY[] + cos*XZ[];                                                 \
    AT[] = cos*XY[] - sin*XZ[];                                                 \
                                                                                \
    RA[] = sin*YX[] + cos*ZX[];                                                 \
    RR[] = sin*sin*YY[] + sin*cos*ZY[]                                          \
         + sin*cos*YZ[] + cos*cos*ZZ[];                                         \
    RT[] = sin*cos*YY[] + cos*cos*ZY[]                                          \
         - sin*sin*YZ[] - sin*cos*ZZ[];                                         \
                                                                                \
    TA[] = cos*YX[] - sin*ZX[];                                                 \
    TR[] = sin*cos*YY[] - sin*sin*ZY[]                                          \
         + cos*cos*YZ[] - sin*cos*ZZ[];                                         \
    TT[] = cos*cos*YY[] - sin*cos*ZY[]                                          \
         - sin*cos*YZ[] + sin*sin*ZZ[];                                         \
}

#define TRANSSYMMTENS(AA,AR,AT,RR,RT,TT,XX,XY,XZ,YY,YZ,ZZ)                      \
                                                                                \
foreach()                                                                       \
{                                                                               \
    double sin = y/sqrt(sq(y) + sq(z));                                         \
    double cos = z/sqrt(sq(y) + sq(z));                                         \
                                                                                \
    AA[] = XX[];                                                                \
    AR[] = sin*XY[] + cos*XZ[];                                                 \
    AT[] = cos*XY[] - sin*XZ[];                                                 \
                                                                                \
    RR[] = sin*sin*YY[] + sin*cos*YZ[]                                          \
         + sin*cos*YZ[] + cos*cos*ZZ[];                                         \
    RT[] = sin*cos*YY[] - sin*sin*YZ[]                                          \
         + cos*cos*YZ[] - sin*cos*ZZ[];                                         \
                                                                                \
    TT[] = cos*cos*YY[] - sin*cos*YZ[]                                          \
         - sin*cos*YZ[] + sin*sin*ZZ[];                                         \
}

// Macro to perform surface averaging

#define CREATEDATA_CYLINDRICAL(NFIELDS)                                         \
                                                                                \
    int Nf = NFIELDS;                                                           \
                                                                                \
    scalar * fieldData = NULL;                                                  \
    double * surfData = malloc(sizeof(double)*Nf*Nr*Na);

#define CREATEDATA_CARTESIAN(NFIELDS)                                           \
                                                                                \
    int Nf = NFIELDS;                                                           \
                                                                                \
    scalar * fieldData = NULL;                                                  \
    double * lineData = malloc(sizeof(double)*Nf*Nr);

#define APPENDTODATA(FIELD)                                                     \
                                                                                \
    fieldData = list_append(fieldData, FIELD);

#define WRITEDATA_CYLINDRICAL(PHASENAME,FIELDNAME,WRITECOORDS)                  \
                                                                                \
    surfAverage(surfData,fieldData,aArr,rArr,tArr);                             \
                                                                                \
    if (dirIndex >= 0 && master)                                                \
    {                                                                           \
        writeSpatialDataAscii                                                   \
        (                                                                       \
            surfCoords,                                                         \
            surfData,                                                           \
            STRINGIFY(PHASENAME) "-" STRINGIFY(FIELDNAME) ".txt",               \
            dirIndex,                                                           \
            2,                                                                  \
            Nf,                                                                 \
            Nr*Na,                                                              \
            WRITECOORDS                                                         \
        );                                                                      \
    }                                                                           \
                                                                                \
    free(surfData);                                                             \
    free(fieldData);

#define WRITEDATA_CARTESIAN(PHASENAME,FIELDNAME,WRITECOORDS)                    \
                                                                                \
    lineAverage(lineData,fieldData,xArr,yArr,zArr);                             \
                                                                                \
    if (dirIndex >= 0 && master)                                                \
    {                                                                           \
        writeSpatialDataAscii                                                   \
        (                                                                       \
            lineCoords,                                                         \
            lineData,                                                           \
            STRINGIFY(PHASENAME) "-" STRINGIFY(FIELDNAME) ".txt",               \
            dirIndex,                                                           \
            1,                                                                  \
            Nf,                                                                 \
            Nr,                                                                 \
            WRITECOORDS                                                         \
        );                                                                      \
    }                                                                           \
                                                                                \
    free(lineData);                                                             \
    free(fieldData);

// Macro to update the average with an old mean weight (wa) and new field weight
// (wb). By setting wa to zero and wb to unity, the mean is reset

#define UPDATEAVERAGE(MEAN,FIELD) MEAN = wa*(MEAN) + wb*(FIELD)

// Macro to create a variable with the phasename prepended

#define GLUE2(PHASENAME,FIELDNAME) PHASENAME##FIELDNAME
#define GLUE(PHASENAME,FIELDNAME) GLUE2(PHASENAME,FIELDNAME)

// Macros to create a string from macro variable

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

// Macro to create mean fields for phase PHASENAME, used for temporal averaging

#define CREATEMEANFIELDS(PHASENAME)                                             \
                                                                                \
scalar PHASENAME##fMean[];                                                      \
                                                                                \
scalar PHASENAME##ffMean[];                                                     \
                                                                                \
vector PHASENAME##fuiMean[];                                                    \
                                                                                \
vector PHASENAME##fuxujMean[];                                                  \
vector PHASENAME##fuyujMean[];                                                  \
vector PHASENAME##fuzujMean[];                                                  \
                                                                                \
vector PHASENAME##fduxdjMean[];                                                 \
vector PHASENAME##fduydjMean[];                                                 \
vector PHASENAME##fduzdjMean[];                                                 \
                                                                                \
scalar PHASENAME##fuiuiMean[];                                                  \
                                                                                \
vector PHASENAME##fuiuiujMean[];                                                \
                                                                                \
scalar PHASENAME##fQMean[];                                                     \
                                                                                \
vector PHASENAME##fQujMean[];                                                   \
                                                                                \
scalar PHASENAME##dfQujdjMean2[];                                               \
                                                                                \
scalar PHASENAME##dfuidjnuduidjMean[];                                          \
                                                                                \
scalar PHASENAME##fuidnuduidjdjMean[];                                          \
                                                                                \
vector PHASENAME##fuinudujdiMean[];                                             \
                                                                                \
vector PHASENAME##fnuuiSijMean[];                                               \
                                                                                \
vector PHASENAME##dfdjnuduidjMean[];                                            \
                                                                                \
vector PHASENAME##fdnuduidjdjMean[];                                            \
                                                                                \
vector PHASENAME##fnuduidxMean[];                                               \
vector PHASENAME##fnuduidyMean[];                                               \
vector PHASENAME##fnuduidzMean[];                                               \
                                                                                \
vector PHASENAME##fnuSxjMean[];                                                 \
vector PHASENAME##fnuSyjMean[];                                                 \
vector PHASENAME##fnuSzjMean[];                                                 \
                                                                                \
scalar PHASENAME##fnuSijduidjMean[];                                            \
                                                                                \
scalar PHASENAME##fTxxMean[];                                                   \
scalar PHASENAME##fTxyMean[];                                                   \
scalar PHASENAME##fTxzMean[];                                                   \
scalar PHASENAME##fTyyMean[];                                                   \
scalar PHASENAME##fTyzMean[];                                                   \
scalar PHASENAME##fTzzMean[];                                                   \
                                                                                \
vector PHASENAME##rhoTijdGdjMean[];                                             \
scalar PHASENAME##rhoTijdGdjuiMean[];                                           \
                                                                                \
vector PHASENAME##fsiMean[];                                                    \
scalar PHASENAME##fuisiMean[];

// Macro to disable storage of mean fields in dump file

#define SETNODUMPMEANFIELDS(PHASENAME)                                          \
                                                                                \
PHASENAME##fMean.nodump = true;                                                 \
                                                                                \
PHASENAME##ffMean.nodump = true;                                                \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fuiMean.x.nodump = true;                                         \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fuxujMean.x.nodump = true;                                       \
foreach_dimension()                                                             \
    PHASENAME##fuyujMean.x.nodump = true;                                       \
foreach_dimension()                                                             \
    PHASENAME##fuzujMean.x.nodump = true;                                       \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fduxdjMean.x.nodump = true;                                      \
foreach_dimension()                                                             \
    PHASENAME##fduydjMean.x.nodump = true;                                      \
foreach_dimension()                                                             \
    PHASENAME##fduzdjMean.x.nodump = true;                                      \
                                                                                \
PHASENAME##fuiuiMean.nodump = true;                                             \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fuiuiujMean.x.nodump = true;                                     \
                                                                                \
PHASENAME##fQMean.nodump = true;                                                \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fQujMean.x.nodump = true;                                        \
                                                                                \
PHASENAME##dfQujdjMean2.nodump = true;                                          \
                                                                                \
PHASENAME##dfuidjnuduidjMean.nodump = true;                                     \
                                                                                \
PHASENAME##fuidnuduidjdjMean.nodump = true;                                     \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fuinudujdiMean.x.nodump = true;                                  \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fnuuiSijMean.x.nodump = true;                                    \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##dfdjnuduidjMean.x.nodump = true;                                 \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fdnuduidjdjMean.x.nodump = true;                                 \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fnuduidxMean.x.nodump = true;                                    \
foreach_dimension()                                                             \
    PHASENAME##fnuduidyMean.x.nodump = true;                                    \
foreach_dimension()                                                             \
    PHASENAME##fnuduidzMean.x.nodump = true;                                    \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fnuSxjMean.x.nodump = true;                                      \
foreach_dimension()                                                             \
    PHASENAME##fnuSyjMean.x.nodump = true;                                      \
foreach_dimension()                                                             \
    PHASENAME##fnuSzjMean.x.nodump = true;                                      \
                                                                                \
PHASENAME##fnuSijduidjMean.nodump = true;                                       \
                                                                                \
PHASENAME##fTxxMean.nodump = true;                                              \
PHASENAME##fTxyMean.nodump = true;                                              \
PHASENAME##fTxzMean.nodump = true;                                              \
PHASENAME##fTyyMean.nodump = true;                                              \
PHASENAME##fTyzMean.nodump = true;                                              \
PHASENAME##fTzzMean.nodump = true;                                              \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##rhoTijdGdjMean.x.nodump = true;                                  \
                                                                                \
PHASENAME##rhoTijdGdjuiMean.nodump = true;                                      \
                                                                                \
foreach_dimension()                                                             \
    PHASENAME##fsiMean.x.nodump = true;                                         \
PHASENAME##fuisiMean.nodump = true;

// Macro to perform the averaging of the given phase

#define PHASEAVERAGE(PHASENAME)                                                 \
{                                                                               \
    /* Compute phase-specific gradients */                                      \
                                                                                \
    scalar rho = new scalar;                                                    \
    scalar Gamma = new scalar;                                                  \
    scalar Q = new scalar;                                                      \
    scalar fQ = new scalar;                                                     \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        rho[] = rho(f[]);                                                       \
        Gamma[] = PHASENAME[]/rho[];                                            \
        Q[] = p[]/rho[] + (G.x*(x-Z.x)+G.y*(y-Z.y)+G.z*(z-Z.z));                \
        fQ[] = PHASENAME[]*Q[];                                                 \
    }                                                                           \
                                                                                \
    vector dGdj = new vector;                                                   \
    vector dfdj = new vector;                                                   \
    vector dfQdj = new vector;                                                  \
                                                                                \
    gradzg(dGdj,Gamma);                                                         \
    gradzg(dfdj,PHASENAME);                                                     \
    gradzg(dfQdj,fQ);                                                           \
                                                                                \
    delete({Gamma, fQ});                                                        \
                                                                                \
    /* Update phase-specific averages */                                        \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        /* Get mixture properties. Add static pressure. */                      \
                                                                                \
        const double nu = mu(f[])/rho[];                                        \
        const double fk = PHASENAME[];                                          \
                                                                                \
        const double duxdx = duxdj.x[];                                         \
        const double duxdy = duxdj.y[];                                         \
        const double duxdz = duxdj.z[];                                         \
        const double duydx = duydj.x[];                                         \
        const double duydy = duydj.y[];                                         \
        const double duydz = duydj.z[];                                         \
        const double duzdx = duzdj.x[];                                         \
        const double duzdy = duzdj.y[];                                         \
        const double duzdz = duzdj.z[];                                         \
                                                                                \
        const double Sxx = 0.5*(duxdx+duxdx);                                   \
        const double Sxy = 0.5*(duxdy+duydx);                                   \
        const double Sxz = 0.5*(duxdz+duzdx);                                   \
        const double Syx = 0.5*(duydx+duxdy);                                   \
        const double Syy = 0.5*(duydy+duydy);                                   \
        const double Syz = 0.5*(duydz+duzdy);                                   \
        const double Szx = 0.5*(duzdx+duxdz);                                   \
        const double Szy = 0.5*(duzdy+duydz);                                   \
        const double Szz = 0.5*(duzdz+duzdz);                                   \
                                                                                \
        const double Txx = 2.0*nu*Sxx - Q[];                                    \
        const double Txy = 2.0*nu*Sxy;                                          \
        const double Txz = 2.0*nu*Sxz;                                          \
        const double Tyx = 2.0*nu*Syx;                                          \
        const double Tyy = 2.0*nu*Syy - Q[];                                    \
        const double Tyz = 2.0*nu*Syz;                                          \
        const double Tzx = 2.0*nu*Szx;                                          \
        const double Tzy = 2.0*nu*Szy;                                          \
        const double Tzz = 2.0*nu*Szz - Q[];                                    \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fMean[], fk);                                  \
        UPDATEAVERAGE(PHASENAME##ffMean[], fk*fk);                              \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fuiMean.x[], fk*u.x[]);                        \
        UPDATEAVERAGE(PHASENAME##fuiMean.y[], fk*u.y[]);                        \
        UPDATEAVERAGE(PHASENAME##fuiMean.z[], fk*u.z[]);                        \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fuxujMean.x[], fk*u.x[]*u.x[]);                \
        UPDATEAVERAGE(PHASENAME##fuxujMean.y[], fk*u.x[]*u.y[]);                \
        UPDATEAVERAGE(PHASENAME##fuxujMean.z[], fk*u.x[]*u.z[]);                \
        UPDATEAVERAGE(PHASENAME##fuyujMean.x[], fk*u.y[]*u.x[]);                \
        UPDATEAVERAGE(PHASENAME##fuyujMean.y[], fk*u.y[]*u.y[]);                \
        UPDATEAVERAGE(PHASENAME##fuyujMean.z[], fk*u.y[]*u.z[]);                \
        UPDATEAVERAGE(PHASENAME##fuzujMean.x[], fk*u.z[]*u.x[]);                \
        UPDATEAVERAGE(PHASENAME##fuzujMean.y[], fk*u.z[]*u.y[]);                \
        UPDATEAVERAGE(PHASENAME##fuzujMean.z[], fk*u.z[]*u.z[]);                \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fduxdjMean.x[], fk*duxdx);                     \
        UPDATEAVERAGE(PHASENAME##fduxdjMean.y[], fk*duxdy);                     \
        UPDATEAVERAGE(PHASENAME##fduxdjMean.z[], fk*duxdz);                     \
        UPDATEAVERAGE(PHASENAME##fduydjMean.x[], fk*duydx);                     \
        UPDATEAVERAGE(PHASENAME##fduydjMean.y[], fk*duydy);                     \
        UPDATEAVERAGE(PHASENAME##fduydjMean.z[], fk*duydz);                     \
        UPDATEAVERAGE(PHASENAME##fduzdjMean.x[], fk*duzdx);                     \
        UPDATEAVERAGE(PHASENAME##fduzdjMean.y[], fk*duzdy);                     \
        UPDATEAVERAGE(PHASENAME##fduzdjMean.z[], fk*duzdz);                     \
                                                                                \
        double uiui = sq(u.x[]) + sq(u.y[]) + sq(u.z[]);                        \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuiuiMean[],                                             \
            fk*uiui                                                             \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuiuiujMean.x[],                                         \
            fk*uiui*u.x[]                                                       \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuiuiujMean.y[],                                         \
            fk*uiui*u.y[]                                                       \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuiuiujMean.z[],                                         \
            fk*uiui*u.z[]                                                       \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fQMean[], fk*Q[]);                             \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fQujMean.x[], fk*Q[]*u.x[]);                   \
        UPDATEAVERAGE(PHASENAME##fQujMean.y[], fk*Q[]*u.y[]);                   \
        UPDATEAVERAGE(PHASENAME##fQujMean.z[], fk*Q[]*u.z[]);                   \
                                                                                \
        delete({Q});                                                            \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##dfQujdjMean2[],                                          \
            u.x[]*dfQdj.x[]                                                     \
          + u.y[]*dfQdj.y[]                                                     \
          + u.z[]*dfQdj.z[]                                                     \
        );                                                                      \
                                                                                \
        delete((scalar*){dfQdj});                                               \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##dfuidjnuduidjMean[],                                     \
            nu                                                                  \
          * (                                                                   \
                (fk*duxdx+u.x[]*dfdj.x[])*duxdx                                 \
              + (fk*duxdy+u.x[]*dfdj.y[])*duxdy                                 \
              + (fk*duxdz+u.x[]*dfdj.z[])*duxdz                                 \
              + (fk*duydx+u.y[]*dfdj.x[])*duydx                                 \
              + (fk*duydy+u.y[]*dfdj.y[])*duydy                                 \
              + (fk*duydz+u.y[]*dfdj.z[])*duydz                                 \
              + (fk*duzdx+u.z[]*dfdj.x[])*duzdx                                 \
              + (fk*duzdy+u.z[]*dfdj.y[])*duzdy                                 \
              + (fk*duzdz+u.z[]*dfdj.z[])*duzdz                                 \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuidnuduidjdjMean[],                                     \
            fk                                                                  \
          * (                                                                   \
                u.x[]*dnuduidjdj.x[]                                            \
              + u.y[]*dnuduidjdj.y[]                                            \
              + u.z[]*dnuduidjdj.z[]                                            \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuinudujdiMean.x[],                                      \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*duxdx                                                     \
              + u.y[]*duxdy                                                     \
              + u.z[]*duxdz                                                     \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuinudujdiMean.y[],                                      \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*duydx                                                     \
              + u.y[]*duydy                                                     \
              + u.z[]*duydz                                                     \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuinudujdiMean.z[],                                      \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*duzdx                                                     \
              + u.y[]*duzdy                                                     \
              + u.z[]*duzdz                                                     \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fnuuiSijMean.x[],                                        \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*Sxx                                                       \
              + u.y[]*Syx                                                       \
              + u.z[]*Szx                                                       \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fnuuiSijMean.y[],                                        \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*Sxy                                                       \
              + u.y[]*Syy                                                       \
              + u.z[]*Szy                                                       \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fnuuiSijMean.z[],                                        \
            fk*nu                                                               \
          * (                                                                   \
                u.x[]*Sxz                                                       \
              + u.y[]*Syz                                                       \
              + u.z[]*Szz                                                       \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##dfdjnuduidjMean.x[],                                     \
            dfdj.x[]*nu*duxdx + dfdj.y[]*nu*duxdy + dfdj.z[]*nu*duxdz           \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##dfdjnuduidjMean.y[],                                     \
            dfdj.x[]*nu*duydx + dfdj.y[]*nu*duydy + dfdj.z[]*nu*duydz           \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##dfdjnuduidjMean.z[],                                     \
            dfdj.x[]*nu*duzdx + dfdj.y[]*nu*duzdy + dfdj.z[]*nu*duzdz           \
        );                                                                      \
                                                                                \
        delete((scalar*){dfdj});                                                \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fdnuduidjdjMean.x[],                                     \
            fk*dnuduidjdj.x[]                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fdnuduidjdjMean.y[],                                     \
            fk*dnuduidjdj.y[]                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fdnuduidjdjMean.z[],                                     \
            fk*dnuduidjdj.z[]                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fnuduidxMean.x[], fk*nu*duxdx);                \
        UPDATEAVERAGE(PHASENAME##fnuduidxMean.y[], fk*nu*duydx);                \
        UPDATEAVERAGE(PHASENAME##fnuduidxMean.z[], fk*nu*duzdx);                \
        UPDATEAVERAGE(PHASENAME##fnuduidyMean.x[], fk*nu*duxdy);                \
        UPDATEAVERAGE(PHASENAME##fnuduidyMean.y[], fk*nu*duydy);                \
        UPDATEAVERAGE(PHASENAME##fnuduidyMean.z[], fk*nu*duzdy);                \
        UPDATEAVERAGE(PHASENAME##fnuduidzMean.x[], fk*nu*duxdz);                \
        UPDATEAVERAGE(PHASENAME##fnuduidzMean.y[], fk*nu*duydz);                \
        UPDATEAVERAGE(PHASENAME##fnuduidzMean.z[], fk*nu*duzdz);                \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fnuSxjMean.x[], fk*nu*Sxx);                    \
        UPDATEAVERAGE(PHASENAME##fnuSxjMean.y[], fk*nu*Sxy);                    \
        UPDATEAVERAGE(PHASENAME##fnuSxjMean.z[], fk*nu*Sxz);                    \
        UPDATEAVERAGE(PHASENAME##fnuSyjMean.x[], fk*nu*Syx);                    \
        UPDATEAVERAGE(PHASENAME##fnuSyjMean.y[], fk*nu*Syy);                    \
        UPDATEAVERAGE(PHASENAME##fnuSyjMean.z[], fk*nu*Syz);                    \
        UPDATEAVERAGE(PHASENAME##fnuSzjMean.x[], fk*nu*Szx);                    \
        UPDATEAVERAGE(PHASENAME##fnuSzjMean.y[], fk*nu*Szy);                    \
        UPDATEAVERAGE(PHASENAME##fnuSzjMean.z[], fk*nu*Szz);                    \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fnuSijduidjMean[],                                       \
            fk*nu                                                               \
          * (                                                                   \
                Sxx*duxdx                                                       \
              + Sxy*duxdy                                                       \
              + Sxz*duxdz                                                       \
              + Syx*duydx                                                       \
              + Syy*duydy                                                       \
              + Syz*duydz                                                       \
              + Szx*duzdx                                                       \
              + Szy*duzdy                                                       \
              + Szz*duzdz                                                       \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fTxxMean[], fk*Txx);                           \
        UPDATEAVERAGE(PHASENAME##fTxyMean[], fk*Txy);                           \
        UPDATEAVERAGE(PHASENAME##fTxzMean[], fk*Txz);                           \
        UPDATEAVERAGE(PHASENAME##fTyyMean[], fk*Tyy);                           \
        UPDATEAVERAGE(PHASENAME##fTyzMean[], fk*Tyz);                           \
        UPDATEAVERAGE(PHASENAME##fTzzMean[], fk*Tzz);                           \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##rhoTijdGdjMean.x[],                                      \
            rho[]                                                               \
          * (                                                                   \
                Txx*dGdj.x[]                                                    \
              + Txy*dGdj.y[]                                                    \
              + Txz*dGdj.z[]                                                    \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##rhoTijdGdjMean.y[],                                      \
            rho[]                                                               \
          * (                                                                   \
                Tyx*dGdj.x[]                                                    \
              + Tyy*dGdj.y[]                                                    \
              + Tyz*dGdj.z[]                                                    \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##rhoTijdGdjMean.z[],                                      \
            rho[]                                                               \
          * (                                                                   \
                Tzx*dGdj.x[]                                                    \
              + Tzy*dGdj.y[]                                                    \
              + Tzz*dGdj.z[]                                                    \
            )                                                                   \
        );                                                                      \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##rhoTijdGdjuiMean[],                                      \
            rho[]                                                               \
          * (                                                                   \
                                                                                \
              + Txx*dGdj.x[]*u.x[]                                              \
              + Txy*dGdj.y[]*u.x[]                                              \
              + Txz*dGdj.z[]*u.x[]                                              \
              + Tyx*dGdj.x[]*u.y[]                                              \
              + Tyy*dGdj.y[]*u.y[]                                              \
              + Tyz*dGdj.z[]*u.y[]                                              \
              + Tzx*dGdj.x[]*u.z[]                                              \
              + Tzy*dGdj.y[]*u.z[]                                              \
              + Tzz*dGdj.z[]*u.z[]                                              \
            )                                                                   \
         );                                                                     \
                                                                                \
        delete((scalar*){dGdj, rho});                                           \
                                                                                \
        UPDATEAVERAGE(PHASENAME##fsiMean.x[], fk*surf.x[]);                     \
        UPDATEAVERAGE(PHASENAME##fsiMean.y[], fk*surf.y[]);                     \
        UPDATEAVERAGE(PHASENAME##fsiMean.z[], fk*surf.z[]);                     \
                                                                                \
        UPDATEAVERAGE                                                           \
        (                                                                       \
            PHASENAME##fuisiMean[],                                             \
            fk*(u.x[]*surf.x[] + u.y[]*surf.y[] + u.z[]*surf.z[])               \
        );                                                                      \
    }                                                                           \
}

// Macro to take a snapshot of the averages of given phase, including spatial
// averaging

#define PHASESNAPSHOT_CYLINDRICAL(PHASENAME)                                    \
{                                                                               \
    /* Create divergence fields */                                              \
                                                                                \
    scalar dfuiuiujdjMean[];                                                    \
    scalar dfQujdjMean1[];                                                      \
    scalar dfnuuiSijdjMean2[];                                                  \
    scalar dfuinudujdidjMean[];                                                 \
                                                                                \
    /* Compute the required gradients and divergences */                        \
                                                                                \
    divergence(dfuiuiujdjMean, PHASENAME##fuiuiujMean);                         \
    divergence(dfQujdjMean1, PHASENAME##fQujMean);                              \
    divergence(dfnuuiSijdjMean2, PHASENAME##fnuuiSijMean);                      \
    divergence(dfuinudujdidjMean, PHASENAME##fuinudujdiMean);                   \
                                                                                \
    /* Create scalar collection fields */                                       \
                                                                                \
    scalar dfnuuiSijdjMean1[];                                                  \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        dfnuuiSijdjMean1[] =                                                    \
            0.5*PHASENAME##dfuidjnuduidjMean[]                                  \
          + 0.5*PHASENAME##fuidnuduidjdjMean[]                                  \
          + 0.5*dfuinudujdidjMean[];                                            \
    }                                                                           \
                                                                                \
    /* Create cylindrical terms */                                              \
                                                                                \
    scalar dfdaMean[];                                                          \
    scalar dfdrMean[];                                                          \
    scalar dfdtMean[];                                                          \
                                                                                \
    scalar fuaMean[];                                                           \
    scalar furMean[];                                                           \
    scalar futMean[];                                                           \
                                                                                \
    scalar fuauaMean[];                                                         \
    scalar fuaurMean[];                                                         \
    scalar fuautMean[];                                                         \
    scalar fururMean[];                                                         \
    scalar furutMean[];                                                         \
    scalar fututMean[];                                                         \
                                                                                \
    scalar dfuadaMean[];                                                        \
    scalar dfurdaMean[];                                                        \
    scalar dfutdaMean[];                                                        \
    scalar dfuadrMean[];                                                        \
    scalar dfurdrMean[];                                                        \
    scalar dfutdrMean[];                                                        \
    scalar dfuadtMean[];                                                        \
    scalar dfurdtMean[];                                                        \
    scalar dfutdtMean[];                                                        \
                                                                                \
    scalar fduadaMean[];                                                        \
    scalar fdurdaMean[];                                                        \
    scalar fdutdaMean[];                                                        \
    scalar fduadrMean[];                                                        \
    scalar fdurdrMean[];                                                        \
    scalar fdutdrMean[];                                                        \
    scalar fduadtMean[];                                                        \
    scalar fdurdtMean[];                                                        \
    scalar fdutdtMean[];                                                        \
                                                                                \
    scalar dfuiuidaMean[];                                                      \
    scalar dfuiuidrMean[];                                                      \
    scalar dfuiuidtMean[];                                                      \
                                                                                \
    scalar dfuaujdjMean[];                                                      \
    scalar dfurujdjMean[];                                                      \
    scalar dfutujdjMean[];                                                      \
                                                                                \
    scalar dfQdaMean[];                                                         \
    scalar dfQdrMean[];                                                         \
    scalar dfQdtMean[];                                                         \
                                                                                \
    scalar dfnuSajdjMean1[];                                                    \
    scalar dfnuSrjdjMean1[];                                                    \
    scalar dfnuStjdjMean1[];                                                    \
                                                                                \
    scalar dfnuSajdjMean2[];                                                    \
    scalar dfnuSrjdjMean2[];                                                    \
    scalar dfnuStjdjMean2[];                                                    \
                                                                                \
    scalar fnuSaaMean[];                                                        \
    scalar fnuSarMean[];                                                        \
    scalar fnuSatMean[];                                                        \
    scalar fnuSrrMean[];                                                        \
    scalar fnuSrtMean[];                                                        \
    scalar fnuSttMean[];                                                        \
                                                                                \
    scalar fTaaMean[];                                                          \
    scalar fTarMean[];                                                          \
    scalar fTatMean[];                                                          \
    scalar fTrrMean[];                                                          \
    scalar fTrtMean[];                                                          \
    scalar fTttMean[];                                                          \
                                                                                \
    scalar rhoTajdGdjMean[];                                                    \
    scalar rhoTrjdGdjMean[];                                                    \
    scalar rhoTtjdGdjMean[];                                                    \
                                                                                \
    scalar fsaMean[];                                                           \
    scalar fsrMean[];                                                           \
    scalar fstMean[];                                                           \
                                                                                \
    /* Compute the cylindrical terms by transformation */                       \
                                                                                \
    vector dfdjMean = new vector;                                               \
                                                                                \
    grad(dfdjMean, PHASENAME##fMean);                                           \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfdaMean,                                                               \
        dfdrMean,                                                               \
        dfdtMean,                                                               \
        dfdjMean.x,                                                             \
        dfdjMean.y,                                                             \
        dfdjMean.z                                                              \
    )                                                                           \
                                                                                \
    delete((scalar*){dfdjMean});                                                \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        fuaMean,                                                                \
        furMean,                                                                \
        futMean,                                                                \
        PHASENAME##fuiMean.x,                                                   \
        PHASENAME##fuiMean.y,                                                   \
        PHASENAME##fuiMean.z                                                    \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSSYMMTENS                                                               \
    (                                                                           \
        fuauaMean,                                                              \
        fuaurMean,                                                              \
        fuautMean,                                                              \
        fururMean,                                                              \
        furutMean,                                                              \
        fututMean,                                                              \
        PHASENAME##fuxujMean.x,                                                 \
        PHASENAME##fuxujMean.y,                                                 \
        PHASENAME##fuxujMean.z,                                                 \
        PHASENAME##fuyujMean.y,                                                 \
        PHASENAME##fuyujMean.z,                                                 \
        PHASENAME##fuzujMean.z                                                  \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfuxdjMean = new vector;                                             \
    vector dfuydjMean = new vector;                                             \
    vector dfuzdjMean = new vector;                                             \
                                                                                \
    grad(dfuxdjMean, PHASENAME##fuiMean.x);                                     \
    grad(dfuydjMean, PHASENAME##fuiMean.y);                                     \
    grad(dfuzdjMean, PHASENAME##fuiMean.z);                                     \
                                                                                \
    TRANSTENS                                                                   \
    (                                                                           \
        dfuadaMean,                                                             \
        dfuadrMean,                                                             \
        dfuadtMean,                                                             \
        dfurdaMean,                                                             \
        dfurdrMean,                                                             \
        dfurdtMean,                                                             \
        dfutdaMean,                                                             \
        dfutdrMean,                                                             \
        dfutdtMean,                                                             \
        dfuxdjMean.x,                                                           \
        dfuxdjMean.y,                                                           \
        dfuxdjMean.z,                                                           \
        dfuydjMean.x,                                                           \
        dfuydjMean.y,                                                           \
        dfuydjMean.z,                                                           \
        dfuzdjMean.x,                                                           \
        dfuzdjMean.y,                                                           \
        dfuzdjMean.z                                                            \
    )                                                                           \
                                                                                \
    delete((scalar*){dfuxdjMean, dfuydjMean, dfuzdjMean});                      \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSTENS                                                                   \
    (                                                                           \
        fduadaMean,                                                             \
        fduadrMean,                                                             \
        fduadtMean,                                                             \
        fdurdaMean,                                                             \
        fdurdrMean,                                                             \
        fdurdtMean,                                                             \
        fdutdaMean,                                                             \
        fdutdrMean,                                                             \
        fdutdtMean,                                                             \
        PHASENAME##fduxdjMean.x,                                                \
        PHASENAME##fduxdjMean.y,                                                \
        PHASENAME##fduxdjMean.z,                                                \
        PHASENAME##fduydjMean.x,                                                \
        PHASENAME##fduydjMean.y,                                                \
        PHASENAME##fduydjMean.z,                                                \
        PHASENAME##fduzdjMean.x,                                                \
        PHASENAME##fduzdjMean.y,                                                \
        PHASENAME##fduzdjMean.z                                                 \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfuiuidjMean = new vector;                                           \
                                                                                \
    gradzg(dfuiuidjMean, PHASENAME##fuiuiMean);                                 \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfuiuidaMean,                                                           \
        dfuiuidrMean,                                                           \
        dfuiuidtMean,                                                           \
        dfuiuidjMean.x,                                                         \
        dfuiuidjMean.y,                                                         \
        dfuiuidjMean.z                                                          \
    )                                                                           \
                                                                                \
    delete((scalar*){dfuiuidjMean});                                            \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfuiujdjMean = new vector;                                           \
                                                                                \
    divergence(dfuiujdjMean.x, PHASENAME##fuxujMean);                           \
    divergence(dfuiujdjMean.y, PHASENAME##fuyujMean);                           \
    divergence(dfuiujdjMean.z, PHASENAME##fuzujMean);                           \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfuaujdjMean,                                                           \
        dfurujdjMean,                                                           \
        dfutujdjMean,                                                           \
        dfuiujdjMean.x,                                                         \
        dfuiujdjMean.y,                                                         \
        dfuiujdjMean.z                                                          \
    )                                                                           \
                                                                                \
    delete((scalar*){dfuiujdjMean});                                            \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfQdjMean = new vector;                                              \
                                                                                \
    gradzg(dfQdjMean, PHASENAME##fQMean);                                       \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfQdaMean,                                                              \
        dfQdrMean,                                                              \
        dfQdtMean,                                                              \
        dfQdjMean.x,                                                            \
        dfQdjMean.y,                                                            \
        dfQdjMean.z                                                             \
    )                                                                           \
                                                                                \
    delete((scalar*){dfQdjMean});                                               \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfnuSijdjMean1 = new vector;                                         \
    vector dfnudujdidjMean = new vector;                                        \
                                                                                \
    divergence(dfnudujdidjMean.x, PHASENAME##fnuduidxMean);                     \
    divergence(dfnudujdidjMean.y, PHASENAME##fnuduidyMean);                     \
    divergence(dfnudujdidjMean.z, PHASENAME##fnuduidzMean);                     \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        foreach_dimension()                                                     \
        {                                                                       \
            dfnuSijdjMean1.x[] =                                                \
                0.5*PHASENAME##dfdjnuduidjMean.x[]                              \
              + 0.5*PHASENAME##fdnuduidjdjMean.x[]                              \
              + 0.5*dfnudujdidjMean.x[];                                        \
        }                                                                       \
    }                                                                           \
                                                                                \
    delete((scalar*){dfnudujdidjMean});                                         \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfnuSajdjMean1,                                                         \
        dfnuSrjdjMean1,                                                         \
        dfnuStjdjMean1,                                                         \
        dfnuSijdjMean1.x,                                                       \
        dfnuSijdjMean1.y,                                                       \
        dfnuSijdjMean1.z                                                        \
    )                                                                           \
                                                                                \
    delete((scalar*){dfnuSijdjMean1});                                          \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfnuSijdjMean2 = new vector;                                         \
                                                                                \
    divergence(dfnuSijdjMean2.x, PHASENAME##fnuSxjMean);                        \
    divergence(dfnuSijdjMean2.y, PHASENAME##fnuSyjMean);                        \
    divergence(dfnuSijdjMean2.z, PHASENAME##fnuSzjMean);                        \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        dfnuSajdjMean2,                                                         \
        dfnuSrjdjMean2,                                                         \
        dfnuStjdjMean2,                                                         \
        dfnuSijdjMean2.x,                                                       \
        dfnuSijdjMean2.y,                                                       \
        dfnuSijdjMean2.z                                                        \
    )                                                                           \
                                                                                \
    delete((scalar*){dfnuSijdjMean2});                                          \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSSYMMTENS                                                               \
    (                                                                           \
        fnuSaaMean,                                                             \
        fnuSarMean,                                                             \
        fnuSatMean,                                                             \
        fnuSrrMean,                                                             \
        fnuSrtMean,                                                             \
        fnuSttMean,                                                             \
        PHASENAME##fnuSxjMean.x,                                                \
        PHASENAME##fnuSxjMean.y,                                                \
        PHASENAME##fnuSxjMean.z,                                                \
        PHASENAME##fnuSyjMean.y,                                                \
        PHASENAME##fnuSyjMean.z,                                                \
        PHASENAME##fnuSzjMean.z                                                 \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSSYMMTENS                                                               \
    (                                                                           \
        fTaaMean,                                                               \
        fTarMean,                                                               \
        fTatMean,                                                               \
        fTrrMean,                                                               \
        fTrtMean,                                                               \
        fTttMean,                                                               \
        PHASENAME##fTxxMean,                                                    \
        PHASENAME##fTxyMean,                                                    \
        PHASENAME##fTxzMean,                                                    \
        PHASENAME##fTyyMean,                                                    \
        PHASENAME##fTyzMean,                                                    \
        PHASENAME##fTzzMean                                                     \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        rhoTajdGdjMean,                                                         \
        rhoTrjdGdjMean,                                                         \
        rhoTtjdGdjMean,                                                         \
        PHASENAME##rhoTijdGdjMean.x,                                            \
        PHASENAME##rhoTijdGdjMean.y,                                            \
        PHASENAME##rhoTijdGdjMean.z                                             \
    )                                                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    TRANSVEC                                                                    \
    (                                                                           \
        fsaMean,                                                                \
        fsrMean,                                                                \
        fstMean,                                                                \
        PHASENAME##fsiMean.x,                                                   \
        PHASENAME##fsiMean.y,                                                   \
        PHASENAME##fsiMean.z                                                    \
    )                                                                           \
                                                                                \
    /* Write 73 fields in three sets, of 31, 28 and 14 fields */                \
    /* The first set contains the coordinates too */                            \
                                                                                \
    {                                                                           \
        CREATEDATA_CYLINDRICAL(31)                                              \
                                                                                \
        APPENDTODATA(PHASENAME##fMean)                                          \
                                                                                \
        APPENDTODATA(dfdaMean)                                                  \
        APPENDTODATA(dfdrMean)                                                  \
        APPENDTODATA(dfdtMean)                                                  \
                                                                                \
        APPENDTODATA(fuaMean)                                                   \
        APPENDTODATA(furMean)                                                   \
        APPENDTODATA(futMean)                                                   \
                                                                                \
        APPENDTODATA(fuauaMean)                                                 \
        APPENDTODATA(fuaurMean)                                                 \
        APPENDTODATA(fuautMean)                                                 \
        APPENDTODATA(fururMean)                                                 \
        APPENDTODATA(furutMean)                                                 \
        APPENDTODATA(fututMean)                                                 \
                                                                                \
        APPENDTODATA(dfuadaMean)                                                \
        APPENDTODATA(dfurdaMean)                                                \
        APPENDTODATA(dfutdaMean)                                                \
        APPENDTODATA(dfuadrMean)                                                \
        APPENDTODATA(dfurdrMean)                                                \
        APPENDTODATA(dfutdrMean)                                                \
        APPENDTODATA(dfuadtMean)                                                \
        APPENDTODATA(dfurdtMean)                                                \
        APPENDTODATA(dfutdtMean)                                                \
                                                                                \
        APPENDTODATA(fduadaMean)                                                \
        APPENDTODATA(fdurdaMean)                                                \
        APPENDTODATA(fdutdaMean)                                                \
        APPENDTODATA(fduadrMean)                                                \
        APPENDTODATA(fdurdrMean)                                                \
        APPENDTODATA(fdutdrMean)                                                \
        APPENDTODATA(fduadtMean)                                                \
        APPENDTODATA(fdurdtMean)                                                \
        APPENDTODATA(fdutdtMean)                                                \
                                                                                \
        WRITEDATA_CYLINDRICAL(PHASENAME,1,true)                                 \
    }                                                                           \
                                                                                \
    {                                                                           \
        CREATEDATA_CYLINDRICAL(28)                                              \
                                                                                \
        APPENDTODATA(dfuiuidaMean)                                              \
        APPENDTODATA(dfuiuidrMean)                                              \
        APPENDTODATA(dfuiuidtMean)                                              \
                                                                                \
        APPENDTODATA(dfuaujdjMean)                                              \
        APPENDTODATA(dfurujdjMean)                                              \
        APPENDTODATA(dfutujdjMean)                                              \
                                                                                \
        APPENDTODATA(dfuiuiujdjMean)                                            \
                                                                                \
        APPENDTODATA(PHASENAME##fQMean)                                         \
                                                                                \
        APPENDTODATA(dfQdaMean)                                                 \
        APPENDTODATA(dfQdrMean)                                                 \
        APPENDTODATA(dfQdtMean)                                                 \
                                                                                \
        APPENDTODATA(dfQujdjMean1)                                              \
        APPENDTODATA(PHASENAME##dfQujdjMean2)                                   \
                                                                                \
        APPENDTODATA(dfnuuiSijdjMean1)                                          \
        APPENDTODATA(dfnuuiSijdjMean2)                                          \
                                                                                \
        APPENDTODATA(dfnuSajdjMean1)                                            \
        APPENDTODATA(dfnuSrjdjMean1)                                            \
        APPENDTODATA(dfnuStjdjMean1)                                            \
                                                                                \
        APPENDTODATA(dfnuSajdjMean2)                                            \
        APPENDTODATA(dfnuSrjdjMean2)                                            \
        APPENDTODATA(dfnuStjdjMean2)                                            \
                                                                                \
        APPENDTODATA(fnuSaaMean)                                                \
        APPENDTODATA(fnuSarMean)                                                \
        APPENDTODATA(fnuSatMean)                                                \
        APPENDTODATA(fnuSrrMean)                                                \
        APPENDTODATA(fnuSrtMean)                                                \
        APPENDTODATA(fnuSttMean)                                                \
                                                                                \
        APPENDTODATA(PHASENAME##fnuSijduidjMean)                                \
                                                                                \
        WRITEDATA_CYLINDRICAL(PHASENAME,2,false)                                \
    }                                                                           \
                                                                                \
    {                                                                           \
        CREATEDATA_CYLINDRICAL(15)                                              \
                                                                                \
        APPENDTODATA(fTaaMean)                                                  \
        APPENDTODATA(fTarMean)                                                  \
        APPENDTODATA(fTatMean)                                                  \
        APPENDTODATA(fTrrMean)                                                  \
        APPENDTODATA(fTrtMean)                                                  \
        APPENDTODATA(fTttMean)                                                  \
                                                                                \
        APPENDTODATA(rhoTajdGdjMean)                                            \
        APPENDTODATA(rhoTrjdGdjMean)                                            \
        APPENDTODATA(rhoTtjdGdjMean)                                            \
                                                                                \
        APPENDTODATA(PHASENAME##rhoTijdGdjuiMean)                               \
                                                                                \
        APPENDTODATA(fsaMean)                                                   \
        APPENDTODATA(fsrMean)                                                   \
        APPENDTODATA(fstMean)                                                   \
                                                                                \
        APPENDTODATA(PHASENAME##fuisiMean)                                      \
                                                                                \
        APPENDTODATA(PHASENAME##ffMean)                                         \
                                                                                \
        WRITEDATA_CYLINDRICAL(PHASENAME,3,false)                                \
    }                                                                           \
                                                                                \
    if (master)                                                                 \
    {                                                                           \
        fprintf                                                                 \
        (                                                                       \
            stdout,                                                             \
            "Wrote averaging data for field %s\n",                              \
            STRINGIFY(PHASENAME)                                                \
        );                                                                      \
                                                                                \
        fflush(stdout);                                                         \
    }                                                                           \
}

#define PHASESNAPSHOT_CARTESIAN(PHASENAME)                                      \
{                                                                               \
    scalar dfuiuiujdjMean[];                                                    \
    scalar dfQujdjMean1[];                                                      \
    scalar dfnuuiSijdjMean2[];                                                  \
    scalar dfuinudujdidjMean[];                                                 \
                                                                                \
    divergence(dfuiuiujdjMean, PHASENAME##fuiuiujMean);                         \
    divergence(dfQujdjMean1, PHASENAME##fQujMean);                              \
    divergence(dfnuuiSijdjMean2, PHASENAME##fnuuiSijMean);                      \
    divergence(dfuinudujdidjMean, PHASENAME##fuinudujdiMean);                   \
                                                                                \
    scalar dfnuuiSijdjMean1[];                                                  \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        dfnuuiSijdjMean1[] =                                                    \
            0.5*PHASENAME##dfuidjnuduidjMean[]                                  \
          + 0.5*PHASENAME##fuidnuduidjdjMean[]                                  \
          + 0.5*dfuinudujdidjMean[];                                            \
    }                                                                           \
                                                                                \
    vector dfuiuidjMean = new vector;                                           \
                                                                                \
    vector dfdjMean = new vector;                                               \
                                                                                \
    grad(dfdjMean, PHASENAME##fMean);                                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    gradzg(dfuiuidjMean, PHASENAME##fuiuiMean);                                 \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfuiujdjMean = new vector;                                           \
                                                                                \
    divergence(dfuiujdjMean.x, PHASENAME##fuxujMean);                           \
    divergence(dfuiujdjMean.y, PHASENAME##fuyujMean);                           \
    divergence(dfuiujdjMean.z, PHASENAME##fuzujMean);                           \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfQdjMean = new vector;                                              \
                                                                                \
    gradzg(dfQdjMean, PHASENAME##fQMean);                                       \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfuxdjMean = new vector;                                             \
    vector dfuydjMean = new vector;                                             \
    vector dfuzdjMean = new vector;                                             \
                                                                                \
    grad(dfuxdjMean, PHASENAME##fuiMean.x);                                     \
    grad(dfuydjMean, PHASENAME##fuiMean.y);                                     \
    grad(dfuzdjMean, PHASENAME##fuiMean.z);                                     \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfnuSijdjMean1 = new vector;                                         \
    vector dfnudujdidjMean = new vector;                                        \
                                                                                \
    divergence(dfnudujdidjMean.x, PHASENAME##fnuduidxMean);                     \
    divergence(dfnudujdidjMean.y, PHASENAME##fnuduidyMean);                     \
    divergence(dfnudujdidjMean.z, PHASENAME##fnuduidzMean);                     \
                                                                                \
    foreach()                                                                   \
    {                                                                           \
        foreach_dimension()                                                     \
        {                                                                       \
            dfnuSijdjMean1.x[] =                                                \
                0.5*PHASENAME##dfdjnuduidjMean.x[]                              \
              + 0.5*PHASENAME##fdnuduidjdjMean.x[]                              \
              + 0.5*dfnudujdidjMean.x[];                                        \
        }                                                                       \
    }                                                                           \
                                                                                \
    delete((scalar*){dfnudujdidjMean});                                         \
                                                                                \
    /* --- */                                                                   \
                                                                                \
    vector dfnuSijdjMean2 = new vector;                                         \
                                                                                \
    divergence(dfnuSijdjMean2.x, PHASENAME##fnuSxjMean);                        \
    divergence(dfnuSijdjMean2.y, PHASENAME##fnuSyjMean);                        \
    divergence(dfnuSijdjMean2.z, PHASENAME##fnuSzjMean);                        \
                                                                                \
    /* Write 73 fields in three sets, of 31, 28 and 14 fields */                \
    /* The first set contains the coordinates too */                            \
                                                                                \
    {                                                                           \
        CREATEDATA_CARTESIAN(31)                                                \
                                                                                \
        APPENDTODATA(PHASENAME##fMean)                                          \
                                                                                \
        APPENDTODATA(dfdjMean.x)                                                \
        APPENDTODATA(dfdjMean.y)                                                \
        APPENDTODATA(dfdjMean.z)                                                \
                                                                                \
        APPENDTODATA(PHASENAME##fuiMean.x)                                      \
        APPENDTODATA(PHASENAME##fuiMean.y)                                      \
        APPENDTODATA(PHASENAME##fuiMean.z)                                      \
                                                                                \
        APPENDTODATA(PHASENAME##fuxujMean.x)                                    \
        APPENDTODATA(PHASENAME##fuxujMean.y)                                    \
        APPENDTODATA(PHASENAME##fuxujMean.z)                                    \
        APPENDTODATA(PHASENAME##fuyujMean.y)                                    \
        APPENDTODATA(PHASENAME##fuyujMean.z)                                    \
        APPENDTODATA(PHASENAME##fuzujMean.z)                                    \
                                                                                \
        APPENDTODATA(dfuxdjMean.x)                                              \
        APPENDTODATA(dfuxdjMean.y)                                              \
        APPENDTODATA(dfuxdjMean.z)                                              \
        APPENDTODATA(dfuydjMean.x)                                              \
        APPENDTODATA(dfuydjMean.y)                                              \
        APPENDTODATA(dfuydjMean.z)                                              \
        APPENDTODATA(dfuzdjMean.x)                                              \
        APPENDTODATA(dfuzdjMean.y)                                              \
        APPENDTODATA(dfuzdjMean.z)                                              \
                                                                                \
        APPENDTODATA(PHASENAME##fduxdjMean.x)                                   \
        APPENDTODATA(PHASENAME##fduxdjMean.y)                                   \
        APPENDTODATA(PHASENAME##fduxdjMean.z)                                   \
        APPENDTODATA(PHASENAME##fduydjMean.x)                                   \
        APPENDTODATA(PHASENAME##fduydjMean.y)                                   \
        APPENDTODATA(PHASENAME##fduydjMean.z)                                   \
        APPENDTODATA(PHASENAME##fduzdjMean.x)                                   \
        APPENDTODATA(PHASENAME##fduzdjMean.y)                                   \
        APPENDTODATA(PHASENAME##fduzdjMean.z)                                   \
                                                                                \
        WRITEDATA_CARTESIAN(PHASENAME,1,true)                                   \
    }                                                                           \
                                                                                \
    {                                                                           \
        CREATEDATA_CARTESIAN(28)                                                \
                                                                                \
        APPENDTODATA(dfuiuidjMean.x)                                            \
        APPENDTODATA(dfuiuidjMean.y)                                            \
        APPENDTODATA(dfuiuidjMean.z)                                            \
                                                                                \
        APPENDTODATA(dfuiujdjMean.x)                                            \
        APPENDTODATA(dfuiujdjMean.y)                                            \
        APPENDTODATA(dfuiujdjMean.z)                                            \
                                                                                \
        APPENDTODATA(dfuiuiujdjMean)                                            \
                                                                                \
        APPENDTODATA(PHASENAME##fQMean)                                         \
                                                                                \
        APPENDTODATA(dfQdjMean.x)                                               \
        APPENDTODATA(dfQdjMean.y)                                               \
        APPENDTODATA(dfQdjMean.z)                                               \
                                                                                \
        APPENDTODATA(dfQujdjMean1)                                              \
        APPENDTODATA(PHASENAME##dfQujdjMean2)                                   \
                                                                                \
        APPENDTODATA(dfnuuiSijdjMean1)                                          \
        APPENDTODATA(dfnuuiSijdjMean2)                                          \
                                                                                \
        APPENDTODATA(dfnuSijdjMean1.x)                                          \
        APPENDTODATA(dfnuSijdjMean1.y)                                          \
        APPENDTODATA(dfnuSijdjMean1.z)                                          \
                                                                                \
        APPENDTODATA(dfnuSijdjMean2.x)                                          \
        APPENDTODATA(dfnuSijdjMean2.y)                                          \
        APPENDTODATA(dfnuSijdjMean2.z)                                          \
                                                                                \
        APPENDTODATA(PHASENAME##fnuSxjMean.x)                                   \
        APPENDTODATA(PHASENAME##fnuSxjMean.y)                                   \
        APPENDTODATA(PHASENAME##fnuSxjMean.z)                                   \
        APPENDTODATA(PHASENAME##fnuSyjMean.y)                                   \
        APPENDTODATA(PHASENAME##fnuSyjMean.z)                                   \
        APPENDTODATA(PHASENAME##fnuSzjMean.z)                                   \
                                                                                \
        APPENDTODATA(PHASENAME##fnuSijduidjMean)                                \
                                                                                \
        WRITEDATA_CARTESIAN(PHASENAME,2,false)                                  \
    }                                                                           \
                                                                                \
    {                                                                           \
        CREATEDATA_CARTESIAN(15)                                                \
                                                                                \
        APPENDTODATA(PHASENAME##fTxxMean)                                       \
        APPENDTODATA(PHASENAME##fTxyMean)                                       \
        APPENDTODATA(PHASENAME##fTxzMean)                                       \
        APPENDTODATA(PHASENAME##fTyyMean)                                       \
        APPENDTODATA(PHASENAME##fTyzMean)                                       \
        APPENDTODATA(PHASENAME##fTzzMean)                                       \
                                                                                \
        APPENDTODATA(PHASENAME##rhoTijdGdjMean.x)                               \
        APPENDTODATA(PHASENAME##rhoTijdGdjMean.y)                               \
        APPENDTODATA(PHASENAME##rhoTijdGdjMean.z)                               \
                                                                                \
        APPENDTODATA(PHASENAME##rhoTijdGdjuiMean)                               \
                                                                                \
        APPENDTODATA(PHASENAME##fsiMean.x)                                      \
        APPENDTODATA(PHASENAME##fsiMean.y)                                      \
        APPENDTODATA(PHASENAME##fsiMean.z)                                      \
                                                                                \
        APPENDTODATA(PHASENAME##fuisiMean)                                      \
                                                                                \
        APPENDTODATA(PHASENAME##ffMean)                                         \
                                                                                \
        WRITEDATA_CARTESIAN(PHASENAME,3,false)                                  \
    }                                                                           \
                                                                                \
    if (master)                                                                 \
    {                                                                           \
        fprintf                                                                 \
        (                                                                       \
            stdout,                                                             \
            "Wrote averaging data for field %s\n",                              \
            STRINGIFY(PHASENAME)                                                \
        );                                                                      \
                                                                                \
        fflush(stdout);                                                         \
    }                                                                           \
}