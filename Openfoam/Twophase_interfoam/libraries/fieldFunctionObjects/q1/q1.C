#include "q1.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q1, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q1,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q1::calc()
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

        const volScalarField sqrMagS(tr(((gradU) & (gradU))));

        const volScalarField q
        (
            0.5*(sqr(tr(gradU)) - sqrMagS)
        );

        const dimensionedScalar smallq(q.dimensions(), 1e-99);

        return store
        (
            resultName_,
            q/stabilise(mag(q)+sqrMagS, smallq)
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

Foam::functionObjects::q1::q1
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "U")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q1::~q1()
{}


// ************************************************************************* //
