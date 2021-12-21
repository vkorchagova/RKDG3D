#include "boundary_integrator_slip.hpp"

// Implementation of class BoundaryIntegratorSlip
BoundaryIntegratorSlip::BoundaryIntegratorSlip(RiemannSolver &rsolver_, const int dim) :
   BoundaryIntegrator(rsolver_,dim) { }

void BoundaryIntegratorSlip::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
   state2 = state1;

   rsolver.Rotate(state2, nor, dim); 

   state2[1] *= -1;

   rsolver.InverseRotate(state2, nor, dim);

};