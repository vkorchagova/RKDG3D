#include "boundary_integrator_open.hpp"
#include "physics.hpp"


// Implementation of class BoundaryIntegratorOpen
BoundaryIntegratorOpen::BoundaryIntegratorOpen(RiemannSolver &rsolver_, const int dim) :
   BoundaryIntegrator(rsolver_,dim) { }

void BoundaryIntegratorOpen::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
   state2 = state1;
};