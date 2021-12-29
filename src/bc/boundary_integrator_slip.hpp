#ifndef BND_INTEGRATOR_SLIP_H
#define BND_INTEGRATOR_SLIP_H

#include "boundary_integrator.hpp"

using namespace std;
using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Slip boundary condition
/// U_n = 0
///
class BoundaryIntegratorSlip : public BoundaryIntegrator
{
private:

public:

    /// Constructor
    BoundaryIntegratorSlip(RiemannSolver &rsolver_, const int dim);

    /// Compute state outside the boundary for Riemann solver
    virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_SLIP_H