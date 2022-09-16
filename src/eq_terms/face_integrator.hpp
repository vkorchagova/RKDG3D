#ifndef FACE_INTEGRATOR_H
#define FACE_INTEGRATOR_H

#include "mfem.hpp"
#include "rs_basic.hpp"


using namespace mfem;

/// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Interior face term: <F.n(u),[w]>
///
class FaceIntegrator : public NonlinearFormIntegrator
{
private:

   /// Riemann solver object
   RiemannSolver& rsolver;

   /// Values of shape function for left FE
   Vector shape1;

   /// Values of shape function for right FE
   Vector shape2;

   /// Values of solution in gaussian points for left FE
   Vector funval1;

   /// Values of solution in gaussian points for right FE
   Vector funval2;

   /// Normal to face
   Vector nor;

   /// Vector to store numerical flux at gaussian point
   Vector fluxN;

   /// Projection of face integration point to the reference space of left FE
   IntegrationPoint eip1;

   /// Projection of face integration point to the reference space of right FE
   IntegrationPoint eip2;

public:

   /// Constructor
   FaceIntegrator(RiemannSolver &rsolver_, const int dim);

   /// Compute part of -<F.n(u), [w]> for the given face 
   virtual void AssembleFaceVector(const FiniteElement &el1,
                                       const FiniteElement &el2,
                                       FaceElementTransformations &Tr,
                                       const Vector &elfun, Vector &elvect);
};

#endif // FACE_INTEGRATOR_H