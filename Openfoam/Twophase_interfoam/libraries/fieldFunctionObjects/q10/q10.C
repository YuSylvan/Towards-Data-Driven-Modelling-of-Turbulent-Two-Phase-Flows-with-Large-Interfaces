#include "q10.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q10, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q10,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q10::calc()
{
    if
    (
        foundObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        )
    )
    {
        const momentumTransportModel& model =
            mesh_.lookupObject<momentumTransportModel>
            (
                momentumTransportModel::typeName
            );

        const volVectorField& U = model.U();
        const surfaceScalarField& phi = model.phi();

        const dimensionedScalar smallU(U.dimensions(), SMALL);

        const volScalarField magDdtU
        (
            mag
            (
                fvc::ddt(U)
              + fvc::div(phi, U)
            )
        );

        // Here DGamma/Ds is approximated as 1/|U|^2*DU/Dt

        const volScalarField q
        (
            magDdtU/sqr(stabilise(mag(U), smallU))
        );

        dimensionedScalar smallq(q.dimensions(), 1E-99);

        return store
        (
            resultName_,
            q/stabilise(mag(q) + 1.0/Lc_, smallq)
        );
    }
    else
    {
        cannotFindObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        );

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::q10::q10
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "U"),
    Lc_("Lc", dimLength, 1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q10::~q10()
{}


// ************************************************************************* //
