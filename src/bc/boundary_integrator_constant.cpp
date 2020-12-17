#include "boundary_integrator_constant.hpp"
#include "physics.hpp"

// Implementation of class BoundaryIntegratorConstant
BoundaryIntegratorConstant::BoundaryIntegratorConstant(RiemannSolver &rsolver_, const int dim, const Vector& _fst) :
   BoundaryIntegrator(rsolver_,dim), fixedState(_fst) { }

void BoundaryIntegratorConstant::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
    state2 = state1;
    rsolver.Rotate(state2, nor, dim); 
    
    for (int i = 0; i < dim+2; ++i)
        state2[i] = 2.0*fixedState[i] - state1[i];

    rsolver.InverseRotate(state2, nor, dim);
};