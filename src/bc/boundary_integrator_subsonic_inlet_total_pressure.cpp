#include "boundary_integrator_subsonic_inlet_total_pressure.hpp"
#include "physics.hpp"

// Implementation of class BoundaryIntegratorSubsonicInletTotalPressure
BoundaryIntegratorSubsonicInletTotalPressure::BoundaryIntegratorSubsonicInletTotalPressure(RiemannSolver &rsolver_, const int dim, double _pTot, double _TTot) :
   BoundaryIntegrator(rsolver_,dim), pTot(_pTot), TTot(_TTot) { }

void BoundaryIntegratorSubsonicInletTotalPressure::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
   state2 = state1;
   //rsolver.Rotate(state2, nor, dim); 

   double rhoU2 = state2[1]*state2[1] + state2[2]*state2[2];
   if (dim == 3)
       rhoU2 += state2[3]*state2[3];

   rhoU2 /= state2[0];

   double MTTotSq = rhoU2 / state2[0] / specific_heat_ratio / gas_constant / TTot;

   double MSq = MTTotSq / (1.0 - 0.5 * (specific_heat_ratio - 1.0) * MTTotSq);

   //totalPressure
   double p = pTot * pow( 1.0 + 0.5 * (specific_heat_ratio - 1.0)*MSq, - specific_heat_ratio / (specific_heat_ratio - 1.0));

   double T = TTot / ( 1.0 + 0.5 * (specific_heat_ratio - 1.0)*MSq );

   double rho = p / gas_constant / T;

   state2[num_equation-1] = p / (specific_heat_ratio - 1.0) + 0.5 * rhoU2 / state2[0] * rho;
   state2[1] = state2[1] / state2[0] * rho;
   state2[2] = state2[2] / state2[0] * rho;
   if (dim == 3)
       state2[3] = state2[3] / state2[0] * rho;
   state2[0] = rho;

   //rsolver.InverseRotate(state2, nor, dim);
};