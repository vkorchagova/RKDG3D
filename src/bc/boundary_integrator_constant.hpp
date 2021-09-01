#ifndef BND_INTEGRATOR_CONSTANT_H
#define BND_INTEGRATOR_CONSTANT_H

#include "boundary_integrator.hpp"

using namespace std;
using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

extern int num_equation;

/// Proc rank 
extern int myRank;

// Interior face term: <F.n(u),[w]>
class BoundaryIntegratorConstant : public BoundaryIntegrator
{
private:

    /// Dirichlet state
    Vector fixedState;

public:

    /// Constructor
    BoundaryIntegratorConstant(RiemannSolver &rsolver_, const int dim, const Vector& _fst);

    // /// Compute part of -<F.n(u), [w]> for the given face 
    // virtual void AssembleFaceVector(const FiniteElement &el1,
    //                                 const FiniteElement &el2,
    //                                 FaceElementTransformations &Tr,
    //                                 const Vector &elfun, Vector &elvect);

    virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_WALL_H