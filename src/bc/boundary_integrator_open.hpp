#ifndef BND_INTEGRATOR_OPEN_H
#define BND_INTEGRATOR_OPEN_H

#include "boundary_integrator.hpp"


using namespace mfem;

/// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Simple open boundary dU/dn = 0
/// Just use equal values inside cell near the boundary and outside
///
class BoundaryIntegratorOpen : public BoundaryIntegrator
{

public:

   /// Constructor
   BoundaryIntegratorOpen(RiemannSolver &rsolver_, const int dim);

   /// Compute state outside the boundary for Riemann solver
   virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_OPEN_H