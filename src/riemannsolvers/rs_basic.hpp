#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "mfem.hpp"


using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Abstract class for Riemann solvers
///
class RiemannSolver
{

protected:

   /// Set of eigenvalues
   Vector lambdaF;

public:
   
   /// Physical flux F, left side
   Vector flux1;

   /// Physical flux F, right side
   Vector flux2;

   /// Primitive variables, left side
   Vector primState1;

   /// Primitive variables, right side
   Vector primState2;

   /// Get normal and tangential components by normal to the cell boundary
   void Rotate(Vector& state, const Vector& nor, int dim);

   /// Get initial components by normal to the cell boundary
   void InverseRotate(Vector& state, const Vector& nor, int dim);

public:

   /// Constructor
   RiemannSolver();

   /// Destructor
   virtual ~RiemannSolver() {};

   /// Compute numerical flux
   virtual double Eval(const Vector &state1, const Vector &state2,
                const Vector &nor, Vector &flux, bool debug = false) = 0;
};

#endif // RIEMANN_SOLVER_H
