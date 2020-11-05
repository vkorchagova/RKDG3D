#ifndef BND_INTEGRATOR_H
#define BND_INTEGRATOR_H

#include "mfem.hpp"
#include "rs_basic.hpp"

using namespace std;
using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

extern const int num_equation;

/// Proc rank 
extern int myRank;

// Interior face term: <F.n(u),[w]>
class BoundaryIntegrator : public NonlinearFormIntegrator
{
protected:

    /// Riemann solver object
    RiemannSolver& rsolver;

    /// Values of shape function for left FE
    Vector shape1;

    /// Values of solution in gaussian points for left FE
    Vector funval1;

    /// Values of "solution" in gaussian points for right FE
    Vector funval2;

    /// Normal to face
    Vector nor;

    /// Vector to store numerical flux at gaussian point
    Vector fluxN;

    /// Projection of face integration point to the reference space of left FE
    IntegrationPoint eip1;

    int dim;

public:

    /// Constructor
    BoundaryIntegrator(RiemannSolver &rsolver_, const int dim);

    /// Compute part of -<F.n(u), [w]> for the given face 
    virtual void AssembleFaceVector(const FiniteElement &el1,
                                    const FiniteElement &el2,
                                    FaceElementTransformations &Tr,
                                    const Vector &elfun, Vector &elvect);

    virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) = 0;
};

#endif // BND_INTEGRATOR_H