#ifndef RIEMANN_SOLVER_HLLC_H
#define RIEMANN_SOLVER_HLLC_H

#include "mfem.hpp"
#include "rs_basic.hpp"


using namespace mfem;

/// Number of equations
extern int num_equation;

///
/// Implements an HLLC flux
///
class RiemannSolverHLLC : public RiemannSolver
{

private:

   /// Aux vector
   Vector D;

   /// F*
   Vector FStar;

   /// U*
   Vector UStar;

   /// Calculate F*
   void getFStar(
       const Vector& state, 
       const Vector& fK, 
       const double rhoL, 
       const double pK, 
       const double SK, 
       const double cK, 
       const double SStar,
       const int dim,
       Vector& FStar
   );

   /// Calculate U*
   void getUStar (
       const Vector& state,
       const double pK,
       const double SK,
       const double cK,
       const double SStar,
       const int dim,
       Vector& UStar
   ) const;

public:

   /// Constructor
   RiemannSolverHLLC();

   /// Destructor
   ~RiemannSolverHLLC() {};

   /// Compute numerical flux
   virtual double Eval(const Vector &state1, const Vector &state2,
                const Vector &nor, Vector &flux, bool debug = false) override;
};

#endif // RIEMANN_SOLVER_HLLC_H
