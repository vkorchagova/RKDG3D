#ifndef BND_INTEGRATOR_SUBSONIC_INLET_TOTAL_PRESSURE_H
#define BND_INTEGRATOR_SUBSONIC_INLET_TOTAL_PRESSURE_H

#include "boundary_integrator.hpp"


using namespace mfem;

// Maximum characteristic speed (updated by integrators)
extern double max_char_speed;

/// Number of equations
extern int num_equation;

extern double specific_heat_ratio;
extern double covolume_constant;
extern double gas_constant;

/// Proc rank 
extern int myRank;

///
/// Subsonic inlet boundary where total pressure and total temperature are fixed
///
class BoundaryIntegratorSubsonicInletTotalPressure : public BoundaryIntegrator
{
private:

   /// Constant value of total pressure
   double pTot;

   /// Constant value of total temperature
   double TTot;

public:

   /// Constructor
   BoundaryIntegratorSubsonicInletTotalPressure(RiemannSolver &rsolver_, const int dim, double _pTot, double _TTot);

   /// Compute state outside the boundary for Riemann solver
   virtual void computeRightState(const Vector& state1, Vector& state2, const Vector& nor) override;
};

#endif // BND_INTEGRATOR_SUBSONIC_INLET_TOTAL_PRESSURE_H