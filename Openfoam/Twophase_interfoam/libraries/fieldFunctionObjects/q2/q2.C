#include "q2.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q2, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q2,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q2::calc()
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

        const volScalarField k(model.k());
        const volVectorField& U = model.U();

        dimensionedScalar smallk(k.dimensions(), 1E-99);

        return store
        (
            resultName_,
            k/stabilise(k+0.5*(U&U), smallk)
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

Foam::functionObjects::q2::q2
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "k")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q2::~q2()
{}


// ************************************************************************* //
