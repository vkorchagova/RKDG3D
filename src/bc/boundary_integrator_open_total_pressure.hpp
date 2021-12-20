#ifndef BND_INTEGRATOR_OPEN_TOTAL_PRESSURE_H
#define BND_INTEGRATOR_OPEN_TOTAL_PRESSURE_H

#include "boundary_integrator.hpp"

using namespace std;
using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

extern int num_equation;

extern double specific_heat_ratio;
extern double covolume_constant;
extern double gas_constant;

/// Proc rank 
extern int myRank;

// Interior face term: <F.n(u),[w]>
class BoundaryIntegratorOpenTotalPressure : public BoundaryIntegrator
{
private:

    double pTotal;
    double TFix;

public:

    /// Constructor
    BoundaryIntegratorOpenTotalPressure(RiemannSolver &rsolver_, const int dim, double _pres = 101325, double _TFix = 293.15);

    // /// Compute part of -<F.n(u), [w]> for the given face 
    // virtual void AssembleFaceVector(const FiniteElement &el1,
    //                                 const FiniteElement &el2,
    //                                 FaceElementTransformations &Tr,
    //                                 const Vector &elfun, Vector &elvect);

    virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_OPEN_H