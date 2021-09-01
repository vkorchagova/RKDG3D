#ifndef RIEMANN_SOLVER_LLF_H
#define RIEMANN_SOLVER_LLF_H

#include "mfem.hpp"
#include "rs_basic.hpp"


using namespace std;
using namespace mfem;

extern int num_equation;

///
/// Implements a simple Rusanov flux
///
class RiemannSolverLLF : public RiemannSolver
{

private:
    Vector lambdaF;

public:

    /// Constructor
    RiemannSolverLLF() : RiemannSolver() {lambdaF.SetSize(num_equation);};

    /// Destructor
    ~RiemannSolverLLF() {};

    /// Compute numerical flux
    virtual double Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, Vector &flux, bool debug = false) override;
};

#endif // RIEMANN_SOLVER_RUS_H
