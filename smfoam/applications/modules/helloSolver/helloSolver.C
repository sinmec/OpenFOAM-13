#include "helloSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "helloLib.H"
#include "IOstreams.H"

#include <cstring>

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(helloSolver, 0);
    addToRunTimeSelectionTable(solver, helloSolver, fvMesh);
}
}

Foam::solvers::helloSolver::helloSolver(fvMesh& mesh)
:
    movingMesh(mesh)
{
    Info<< "[helloSolver] created" << endl;
}


Foam::solvers::helloSolver::~helloSolver()
{}


void Foam::solvers::helloSolver::traceStage(const char* stage) const
{
    Info<< "[helloSolver] " << stage
        << " | time=" << runTime.userTimeName() << endl;
}


void Foam::solvers::helloSolver::preSolve()
{
    traceStage("preSolve");
    movingMesh::preSolve();
}


void Foam::solvers::helloSolver::moveMesh()
{
    traceStage("moveMesh");
    movingMesh::moveMesh();
}


void Foam::solvers::helloSolver::motionCorrector()
{
    traceStage("motionCorrector");
}


void Foam::solvers::helloSolver::prePredictor()
{
    traceStage("prePredictor");
}


void Foam::solvers::helloSolver::momentumTransportPredictor()
{
    traceStage("momentumTransportPredictor");
}


void Foam::solvers::helloSolver::thermophysicalTransportPredictor()
{
    traceStage("thermophysicalTransportPredictor");
}


void Foam::solvers::helloSolver::momentumPredictor()
{
    traceStage("momentumPredictor");
    calculateMolecularHello();
}


void Foam::solvers::helloSolver::thermophysicalPredictor()
{
    traceStage("thermophysicalPredictor");
}


void Foam::solvers::helloSolver::pressureCorrector()
{
    traceStage("pressureCorrector");
}


void Foam::solvers::helloSolver::momentumTransportCorrector()
{
    traceStage("momentumTransportCorrector");
}


void Foam::solvers::helloSolver::thermophysicalTransportCorrector()
{
    traceStage("thermophysicalTransportCorrector");
}


void Foam::solvers::helloSolver::postSolve()
{
    traceStage("postSolve");
}
