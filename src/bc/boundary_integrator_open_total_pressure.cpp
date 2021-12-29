#include "boundary_integrator_open_total_pressure.hpp"
#include "physics.hpp"


// Implementation of class BoundaryIntegratorOpenTotalPressure
BoundaryIntegratorOpenTotalPressure::BoundaryIntegratorOpenTotalPressure(RiemannSolver &rsolver_, const int dim, double _pres, double _TFix) :
   BoundaryIntegrator(rsolver_,dim), pTotal(_pres), TFix(_TFix) { }

void BoundaryIntegratorOpenTotalPressure::computeRightState(const Vector& state1, Vector& state2, const Vector& nor) 
{
    state2 = state1;
    // rsolver.Rotate(state2, nor, dim);

    // double rho2U2 = state2[1]*state2[1] + state2[2]*state2[2];
    // if (dim == 3)
    //     rho2U2 += state2[3]*state2[3];

    // double rhoN = state2[1]*nor[1] + state2[2]*nor[2];
    // if (dim == 3)
    //     rhoN += state2[3]*nor[3];
      
    // double U2 = rho2U2 / state2[0] / state2[0];
    // double epsIn = rhoN > 0 ? (state2[num_equation-1] - 0.5*rho2U2) / state2[0] : gas_constant / (specific_heat_ratio - 1.0) * TFix;

    // double M = ComputeM(state2, dim);

    // double pressure = M > 1 ? pTotal * pow( 1.0 + 0.5 * (specific_heat_ratio - 1.0)*M*M, - specific_heat_ratio / (specific_heat_ratio - 1.0)) : pTotal;
    // double rho = M > 1 ? pressure / (specific_heat_ratio - 1.0) / epsIn : pressure / ( 0.5*U2 + (specific_heat_ratio - 1.0) * epsIn);
    // pressure = M > 1 ? pressure : pressure - 0.5*U2*rho;   

    // double rhoEOut = pressure  / (specific_heat_ratio - 1.0) + 0.5 * rho * U2;

    // state2[num_equation-1] = rhoEOut;
    // state2[1] *= rho / state2[0] ;
    // state2[2] *= rho / state2[0];
    // if (dim == 3)
    //     state2[3] *= rho / state2[0];
    // state2[0] = rho;

    // rsolver.InverseRotate(state2, nor, dim);

};