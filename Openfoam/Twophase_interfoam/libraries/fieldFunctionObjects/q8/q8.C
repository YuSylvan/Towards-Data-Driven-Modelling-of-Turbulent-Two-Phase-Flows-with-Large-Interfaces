#include "q8.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q8, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q8,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q8::calc()
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
        const volScalarField k(model.k());

        const tmp<volTensorField> tgradU(fvc::grad(U));
        const volTensorField& gradU = tgradU();

        const tmp<volVectorField> tgradk(fvc::grad(k));
        const volVectorField& gradk = tgradk();

        const volTensorField S(0.5*(gradU+gradU.T()));

        const volScalarField q(U & gradk);

        const volScalarField norm(mag(model.sigma() && S));

        const dimensionedScalar smallq(q.dimensions(), 1E-99);

        return store
        (
            resultName_,
            q/stabilise(mag(q) + norm, smallq)
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

Foam::functionObjects::q8::q8
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "U")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q8::~q8()
{}


// ************************************************************************* //
