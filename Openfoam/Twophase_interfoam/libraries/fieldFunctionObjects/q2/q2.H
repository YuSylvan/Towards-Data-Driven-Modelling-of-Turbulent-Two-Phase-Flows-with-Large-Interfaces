#ifndef functionObjects_q2_H
#define functionObjects_q2_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q2 Declaration
\*---------------------------------------------------------------------------*/

class q2
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q2");


    // Constructors

        //- Construct from Time and dictionary
        q2
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q2();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
