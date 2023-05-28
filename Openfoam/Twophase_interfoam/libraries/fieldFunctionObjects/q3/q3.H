#ifndef functionObjects_q3_H
#define functionObjects_q3_H

#include "volFields.H"
#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q1 Declaration
\*---------------------------------------------------------------------------*/

class q3
:
    public fieldExpression
{
    // Private Data

        //- Wall distance
        const volScalarField& y_;

    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q3");


    // Constructors

        //- Construct from Time and dictionary
        q3
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q3();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
