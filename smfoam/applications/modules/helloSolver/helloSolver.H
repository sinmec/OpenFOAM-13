#ifndef helloSolver_H
#define helloSolver_H

#include "movingMesh.H"
#include "string.H"

namespace Foam
{
namespace solvers
{

class helloSolver
:
    public movingMesh
{
    void traceStage(const char* stage) const;

public:
    //- Runtime type information
    TypeName("helloSolver");

    // Constructors
    helloSolver(fvMesh& mesh);

    //- Disallow default bitwise copy construction
    helloSolver(const helloSolver&) = delete;

    //- Destructor
    virtual ~helloSolver();

    // Member Functions
    // Estas sao todas as etapas providenciadas pelo foamRun
    // para os solvers.
    virtual void preSolve();
    virtual void moveMesh();
    virtual void motionCorrector();
    virtual void prePredictor();
    virtual void momentumTransportPredictor();
    virtual void thermophysicalTransportPredictor();
    virtual void momentumPredictor();
    virtual void thermophysicalPredictor();
    virtual void pressureCorrector();
    virtual void momentumTransportCorrector();
    virtual void thermophysicalTransportCorrector();
    virtual void postSolve();

    // Member Operators
    //- Disallow default bitwise assignment
    void operator=(const helloSolver&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif