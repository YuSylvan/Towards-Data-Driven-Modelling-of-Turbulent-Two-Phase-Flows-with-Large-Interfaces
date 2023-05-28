#ifndef functionObjects_q6_H
#define functionObjects_q6_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q2 Declaration
\*---------------------------------------------------------------------------*/

class q6
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q6");


    // Constructors

        //- Construct from Time and dictionary
        q6
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q6();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
