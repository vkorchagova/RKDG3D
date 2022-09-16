#include "boundary_integrator_open_fixed_pressure.hpp"
#include "physics.hpp"


// Implementation of class BoundaryIntegratorOpenFixedPressure
BoundaryIntegratorOpenFixedPressure::BoundaryIntegratorOpenFixedPressure(RiemannSolver &rsolver_, const int dim, double _pres) :
   BoundaryIntegrator(rsolver_,dim), pFix(_pres) { }

void BoundaryIntegratorOpenFixedPressure::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
   state2 = state1;
   //rsolver.Rotate(state2, nor, dim);

   double rhoU2 = state2[1]*state2[1] + state2[2]*state2[2];
   if (dim == 3)
       rhoU2 += state2[3]*state2[3];

   double U2 = rhoU2 / state2[0] / state2[0];
   double epsIn = (state2[num_equation-1] - 0.5*rhoU2) / state2[0];

   double rho = pFix / (specific_heat_ratio - 1.0) / epsIn;
   double rhoEOut = pFix  / (specific_heat_ratio - 1.0) + 0.5 * rho * U2;

   state2[num_equation-1] = rhoEOut;
   state2[1] = state1[1] / state1[0] * rho;
   state2[2] = state1[2] / state1[0] * rho;
   if (dim == 3)
       state2[3] = state2[3] / state1[0] * rho;
   state2[0] = rho;

   //rsolver.InverseRotate(state2, nor, dim);
};