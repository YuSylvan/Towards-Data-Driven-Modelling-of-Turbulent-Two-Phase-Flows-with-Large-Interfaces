#include "q6.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(q6, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        q6,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::q6::calc()
{
    if
    (
        foundObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        )
     && foundObject<volScalarField>(fieldName_)
    )
    {
        const momentumTransportModel& model =
            mesh_.lookupObject<momentumTransportModel>
            (
                momentumTransportModel::typeName
            );

        const volScalarField& p = lookupObject<volScalarField>(fieldName_);
        const volVectorField& U = model.U();

        const tmp<volVectorField> tgradP(fvc::grad(p));
        const volVectorField& gradP = tgradP();

        const volVectorField UCmptSqr(cmptMultiply(U,U));

        const tmp<volTensorField> tgradUCmptSqr(fvc::grad(UCmptSqr));
        const volTensorField& gradUCmptSqr = tgradUCmptSqr();

        const volScalarField q(sqrt(gradP & gradP));

        const dimensionedScalar smallq(q.dimensions(), 1E-99);

        const volScalarField norm(0.5*mag(tr(gradUCmptSqr)));

        return store
        (
            resultName_,
            q/stabilise(mag(q) + norm, smallq)
        );
    }
    else
    {
        if (!foundObject<volScalarField>(fieldName_))
        {
            cannotFindObject<volScalarField>(fieldName_);
        }
        else
        {
            cannotFindObject<momentumTransportModel>
            (
                momentumTransportModel::typeName
            );
        }

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::q6::q6
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:

    fieldExpression(name,runTime,dict,typeName,"p")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::q6::~q6()
{}


// ************************************************************************* //
