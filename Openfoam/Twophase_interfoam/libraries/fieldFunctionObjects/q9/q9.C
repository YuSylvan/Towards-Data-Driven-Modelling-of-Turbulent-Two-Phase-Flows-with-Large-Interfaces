#include "q9.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q9, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q9,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q9::calc()
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

        const volScalarField q(mag(model.sigma()));

        const dimensionedScalar smallq(q.dimensions(), 1E-99);

        return store
        (
            resultName_,
            q/stabilise(mag(q) + model.k(), smallq)
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

Foam::functionObjects::q9::q9
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "U")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q9::~q9()
{}


// ************************************************************************* //
