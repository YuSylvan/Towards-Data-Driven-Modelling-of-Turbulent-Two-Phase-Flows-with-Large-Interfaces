// Create a cylindrical grid defined by radial point array rArr, tangential
// point array tArr and axial point array aArr. The radial points are stretched.
// The macro requires a DIAMETER variable which is the diameter of the cylinder
// in unstretched Basilisk coordinates

double rArr[Nr];
double tArr[Nt];
double aArr[Na];

const double rSample = 0.5-0.5*sqrt(2.0)*L0/Na;
const double stretch = 2.0;
const double base = rSample/(stretch*Nr);
const double gf = gradingFactor(base, rSample+base/2.0, Nr);

rArr[0] = rSample;

for (int i = 1; i < Nr; i++)
{
    rArr[i] = rArr[i-1] - (base * pow(gf,i-1) + base * pow(gf,i) )/2.0;
}

for (int i = 0; i < Nt; i++)
{
    tArr[i] = 2.0*pi/Nt*i;
}

for (int i = 0; i < Na; i++)
{
    aArr[i] = L0/Na*(i+0.5);
}

double rCoord[Nr*Na];
double aCoord[Nr*Na];

for (int j = 0; j < Nr; j++)
{
    for (int k = 0; k < Na; k++)
    {
        rCoord[j*Na+k] = rArr[j];
        aCoord[j*Na+k] = aArr[k];
    }
}

double * surfCoords[] = {rCoord, aCoord};
