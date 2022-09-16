#ifndef BND_INTEGRATOR_OPEN_TOTAL_PRESSURE_H
#define BND_INTEGRATOR_OPEN_TOTAL_PRESSURE_H

#include "boundary_integrator.hpp"


using namespace mfem;

/// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

extern double specific_heat_ratio;
extern double covolume_constant;
extern double gas_constant;

/// Proc rank 
extern int myRank;

///
/// Open boundary where total pressure and static temperature are fixed
/// Have a switch between subsonic and supersonic flows
///
class BoundaryIntegratorOpenTotalPressure : public BoundaryIntegrator
{
private:

   /// Constant value of total pressure
   double pTotal;

   /// Constant value of static temperature
   double TFix;

public:

   /// Constructor
   BoundaryIntegratorOpenTotalPressure(RiemannSolver &rsolver_, const int dim, double _pres = 101325, double _TFix = 293.15);

   /// Compute state outside the boundary for Riemann solver
   virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_OPEN_H