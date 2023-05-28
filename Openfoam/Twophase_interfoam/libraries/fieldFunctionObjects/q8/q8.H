#ifndef functionObjects_q8_H
#define functionObjects_q8_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q1 Declaration
\*---------------------------------------------------------------------------*/

class q8
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q8");


    // Constructors

        //- Construct from Time and dictionary
        q8
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q8();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
