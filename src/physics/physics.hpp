#include "mfem.hpp"


using namespace mfem;

/// Physics parameters (updated by case)
extern double specific_heat_ratio;
extern double covolume_constant;
extern double gas_constant;

/// Proc rank 
extern int myRank;

/// Square
template<class T>
inline T sqr(T x) {return x*x;}

// Physicality check (at end)
// Check that the state is physical - enabled in debug mode
bool StateIsPhysical(const Vector &state, const int dim);
bool StateIsPhysicalSay(const Vector &state, const Vector &primState, const int dim);

// Pressure (EOS) computation
double ComputePressure(const Vector &state, int dim);
double ComputeSoundSpeed(double rho, double p);

// Temperature computation
double ComputeTemperature(const Vector &state, int dim);

// Energy by given pressure (EOS)
double ComputeEnergy(double rho, double u, double v, double w, double p);

// Sound speed (EOS) computation
double ComputeSoundSpeed(const Vector &state, int dim);

// Compute the vector flux F(u)
void ComputeFlux(const Vector &state, int dim, DenseMatrix &flux);

// Compute the scalar F(u).n
void ComputeFluxDotN(const Vector &state, const Vector &nor,
                      Vector &fluxN);

void ComputeFluxF(const Vector &state, const Vector &primState, const int dim,
                      Vector &flux);

// Compute the maximum characteristic speed.
double ComputeMaxCharSpeed(const Vector &state, const int dim);

// Compute Mach number
double ComputeM(const Vector &state, const int dim);


// Compute char speeds for defined state
// void ComputeCharSpeedsByState(const Vector &state, Vector& lambdaF, const int dim);

// Compute Einfeldt averaged char speeds via two states
void ComputeEinfeldtCharSpeeds(const Vector &state1, const Vector &state2, const Vector &primState1, const Vector &primState2, Vector& lambdaF, const int dim);

// Compute Toro averaged char speeds via two states
void ComputeToroCharSpeeds(const Vector &state1, const Vector &state2, const Vector &primState1, const Vector &primState2, Vector& lambdaF, const int dim);


// Compute primitive variables (rho, u, v, w, p)
void GetPrimitiveFromConservative(const Vector& state, Vector& primState);

// void ComputeU(const GridFunction& sol, const FiniteElementSpace& vfes, GridFunction& U);
// void ComputeP(const GridFunction& sol, GridFunction& p);