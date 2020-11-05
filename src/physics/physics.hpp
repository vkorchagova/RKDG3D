#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern double specific_heat_ratio;
extern double covolume_constant;

/// Square
template<class T>
inline T sqr(T x) {return x*x;}

// Physicality check (at end)
// Check that the state is physical - enabled in debug mode
bool StateIsPhysical(const Vector &state, const int dim);
bool StateIsPhysicalSay(const Vector &state, const int dim);

// Pressure (EOS) computation
double ComputePressure(const Vector &state, int dim);

// Sound speed (EOS) computation
double ComputeSoundSpeed(const Vector &state, int dim);

// Compute the vector flux F(u)
void ComputeFlux(const Vector &state, int dim, DenseMatrix &flux);

// Compute the scalar F(u).n
void ComputeFluxDotN(const Vector &state, const Vector &nor,
                     Vector &fluxN);

void ComputeFluxF(const Vector &state, const int dim,
                     Vector &flux);

// Compute the maximum characteristic speed.
double ComputeMaxCharSpeed(const Vector &state, const int dim);

// Compute Einfeldt averaged char speeds via two states
void ComputeEinfeldtCharSpeeds(const Vector &state1, const Vector &state2, Vector& lambdaF, const int dim);

// Compute Toro averaged char speeds via two states
void ComputeToroCharSpeeds(const Vector &state1, const Vector &state2, Vector& lambdaF, const int dim);


// Compute primitive variables (rho, u, v, w, p)
void TransformConservativeToPrimitive(Vector& state);