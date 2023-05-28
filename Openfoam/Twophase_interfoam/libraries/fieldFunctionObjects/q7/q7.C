#include "q7.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q7, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q7,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q7::calc()
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

        const tmp<volTensorField> tgradU(fvc::grad(U));
        const volTensorField& gradU = tgradU();

        const volScalarField q(mag((U*U) && gradU));

        const volScalarField norm
        (
            sqrt((U & U) * ((U & gradU) & (U & gradU)))
        );

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

Foam::functionObjects::q7::q7
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "U")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q7::~q7()
{}


// ************************************************************************* //
