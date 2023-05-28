#ifndef functionObjects_q10_H
#define functionObjects_q10_H

#include "fieldExpression.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q1 Declaration
\*---------------------------------------------------------------------------*/

class q10
:
    public fieldExpression
{
    // Private Data

        //- Characteristic length scale of the flow [m]
        dimensionedScalar Lc_;

    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q10");


    // Constructors

        //- Construct from Time and dictionary
        q10
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q10();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
