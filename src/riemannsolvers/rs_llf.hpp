#ifndef RIEMANN_SOLVER_LLF_H
#define RIEMANN_SOLVER_LLF_H

#include "mfem.hpp"
#include "rs_basic.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

///
/// Implements a simple Local Lax -- Fridriechs flux
///
class RiemannSolverLLF : public RiemannSolver
{

public:

    /// Constructor
    RiemannSolverLLF() : RiemannSolver() {lambdaF.SetSize(num_equation);};

    /// Destructor
    ~RiemannSolverLLF() {};

    /// Compute numerical flux
    virtual double Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, Vector &flux, bool debug = false) override;
};

#endif // RIEMANN_SOLVER_LLF_H
