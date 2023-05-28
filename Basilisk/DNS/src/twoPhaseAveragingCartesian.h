// Create a Cartesian grid defined by streamwise point array xArr, wall-normal
// point array yArr and spanwise point array zArr. A LENGTH, HEIGHT and WIDTH
// variable is required. Number of points in x, y and z are set by Na, Nr and
// Nt, for compatibility with the cylindrical averaging.

double xArr[Na];
double yArr[Nr];
double zArr[Nt];

for (int i = 0; i < Na; i++)
{
    xArr[i] = X0 + LENGTH/Na*(i+0.5);
}

for (int i = 0; i < Nr; i++)
{
    yArr[i] = Y0 + HEIGHT/Nr*(i+0.5);
}

for (int i = 0; i < Nt; i++)
{
    zArr[i] = Z0 + WIDTH/Nt*(i+0.5);
}

double * lineCoords[] = {yArr};
