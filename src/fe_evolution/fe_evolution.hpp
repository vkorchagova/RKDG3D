#ifndef FE_EVOLUTION_H
#define FE_EVOLUTION_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

// Problem definition
extern int problem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

extern int num_equation;
extern double specific_heat_ratio;
extern double gas_constant;
extern double covolume_constant;

/// Proc rank 
extern int myRank;

///
/// Time-dependent operator for the right-hand side of the ODE representing the
/// DG weak form
///
class FE_Evolution : public TimeDependentOperator
{
private:

   /// Space dimension (1D, 2D, 3D)
   const int dim;

   /// Finite element space of all variables
   FiniteElementSpace &vfes;

   /// Nonlinear form with boundary and face integrators
   Operator &A;

   /// Matrix of the nonlinear form with domain integrators
   SparseMatrix *Aflux;

   /// Inverse mass matrices
   std::vector<DenseMatrix> Me_inv;

   /// Vector to store solution at point
   mutable Vector state;

   /// Matrix to store physical flux
   mutable DenseMatrix f;

   /// Tensor to store physical flux
   mutable DenseTensor flux;

   /// Vector to store the face terms -<F.n(u), [w]>
   mutable Vector z;

   /// Buffer function to compute physical flux F in gaussian points 
   //  @param state solution at gaussian points
   //  @param flux result
   void GetFlux(const DenseMatrix &state, DenseTensor &flux) const;

public:

   /// Constructor
   FE_Evolution(FiniteElementSpace &_vfes,
                Operator &_A, SparseMatrix *_Aflux);

   /// Compute rhs
   //  @param x solution
   //  @param y result
   virtual void Mult(const Vector &x, Vector &y) const;

   void UpdateInverseMassMatrix();

   /// Destructor
   virtual ~FE_Evolution() { }

   void UpdateAfluxPointer(SparseMatrix *_Aflux) {Aflux = _Aflux; };
};

#endif // FE_EVOLUTION_H