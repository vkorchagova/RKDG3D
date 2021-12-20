#include "boundary_integrator_subsonic_inlet_fixed_pressure.hpp"
#include "physics.hpp"

// Implementation of class BoundaryIntegratorSubsonicInletFixedPressure
BoundaryIntegratorSubsonicInletFixedPressure::BoundaryIntegratorSubsonicInletFixedPressure(RiemannSolver &rsolver_, const int dim, double _pFix, double _TTot) :
   BoundaryIntegrator(rsolver_,dim), pFix(_pFix), TTot(_TTot) { }

void BoundaryIntegratorSubsonicInletFixedPressure::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
    state2 = state1;
    //rsolver.Rotate(state2, nor, dim); 

    double rhoU2 = state2[1]*state2[1] + state2[2]*state2[2];
    if (dim == 3)
        rhoU2 += state2[3]*state2[3];

    rhoU2 /= state2[0];

    double MTTotSq = rhoU2 / state2[0] / specific_heat_ratio / gas_constant / TTot;

    double MSq = MTTotSq / (1.0 - 0.5 * (specific_heat_ratio - 1.0) * MTTotSq);

    // fixedValue
    double p = pFix;

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