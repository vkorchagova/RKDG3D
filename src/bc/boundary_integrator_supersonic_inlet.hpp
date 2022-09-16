#ifndef BND_INTEGRATOR_SUPERSONIC_INLET_H
#define BND_INTEGRATOR_SUPERSONIC_INLET_H

#include "boundary_integrator.hpp"


using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Supersonic inlet. All variables are fixed.
///
class BoundaryIntegratorSupersonicInlet : public BoundaryIntegrator
{
private:

   /// Dirichlet state
   Vector fixedState;

public:

   /// Constructor
   BoundaryIntegratorSupersonicInlet(RiemannSolver &rsolver_, const int dim, const Vector& _fst);

   /// Compute state outside the boundary for Riemann solver
   virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_WALL_H