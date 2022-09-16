#ifndef BND_INTEGRATOR_CHAR_OUTLET_H
#define BND_INTEGRATOR_CHAR_OUTLET_H

#include "boundary_integrator.hpp"


using namespace mfem;

/// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Compute Riemann invariants on boundary
/// !!! TODO direrct computation without Riemann solver !!!
///
class BoundaryIntegratorCharOutlet : public BoundaryIntegrator
{
private:
   /// Dirichlet state
   Vector fixedState;

   double CSpeedOut;
   double MOut;
   double UOut;

public:

   /// Constructor
   BoundaryIntegratorCharOutlet(RiemannSolver &rsolver_, const int dim, const Vector& _fst);

   /// Compute part of -<F.n(u), [w]> for the given face 
   virtual void AssembleFaceVector(const FiniteElement &el1,
                                       const FiniteElement &el2,
                                       FaceElementTransformations &Tr,
                                       const Vector &elfun, Vector &elvect);
   /// Compute state outside the boundary for Riemann solver
   virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_OPEN_H