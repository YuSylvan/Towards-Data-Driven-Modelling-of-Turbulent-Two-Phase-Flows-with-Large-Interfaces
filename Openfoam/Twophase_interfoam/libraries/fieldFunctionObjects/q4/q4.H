#ifndef functionObjects_q4_H
#define functionObjects_q4_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class q4 Declaration
\*---------------------------------------------------------------------------*/

class q4
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the feature field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("q4");


    // Constructors

        //- Construct from Time and dictionary
        q4
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~q4();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
