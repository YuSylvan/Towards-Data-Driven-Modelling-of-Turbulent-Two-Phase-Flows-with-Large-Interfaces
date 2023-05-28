#include "q3.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include "singlePhaseTransportModel.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q3, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q3,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q3::calc()
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

        const singlePhaseTransportModel& transport =
            lookupObject<singlePhaseTransportModel>("transportProperties");

        const volScalarField k(model.k());

        return store
        (
            resultName_,
            min(sqrt(k)*y_/(transport.nu()*50),2.0)
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

Foam::functionObjects::q3::q3
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "k"),
    y_(wallDist::New(this->mesh_).y())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q3::~q3()
{}


// ************************************************************************* //
