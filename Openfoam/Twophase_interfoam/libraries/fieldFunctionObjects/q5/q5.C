#include "q5.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q5, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q5,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q5::calc()
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
        const volScalarField epsilon(model.epsilon());

        const volVectorField& U = model.U();

        const tmp<volTensorField> tgradU(fvc::grad(U));
        const volTensorField& gradU = tgradU();

        const dimensionedScalar smallEps(epsilon.dimensions(), 1E-99);

        const volScalarField q(k/max(epsilon, smallEps));

        const volScalarField magS(0.5*mag(gradU + gradU.T()));

        const dimensionedScalar smallqMagS
        (
            magS.dimensions()*q.dimensions(),
            1E-99
        );

        return store
        (
            resultName_,
            q*magS/stabilise(mag(q)*magS + 1.0, smallqMagS)
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

Foam::functionObjects::q5::q5
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "k")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q5::~q5()
{}


// ************************************************************************* //
