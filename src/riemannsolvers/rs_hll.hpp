#ifndef RIEMANN_SOLVER_HLL_H
#define RIEMANN_SOLVER_HLL_H

#include "mfem.hpp"
#include "rs_basic.hpp"


using namespace mfem;

/// Number of equations
extern int num_equation;

///
/// Implements an HLL flux
///
class RiemannSolverHLL : public RiemannSolver
{

public:

   /// Constructor
   RiemannSolverHLL() : RiemannSolver() {};

   /// Destructor
   ~RiemannSolverHLL() {};

   /// Compute numerical flux
   virtual double Eval(const Vector &state1, const Vector &state2,
                const Vector &nor, Vector &flux, bool debug = false) override;
};

#endif // RIEMANN_SOLVER_HLL_H
