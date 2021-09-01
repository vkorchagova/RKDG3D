#ifndef RIEMANN_SOLVER_RUS_H
#define RIEMANN_SOLVER_RUS_H

#include "mfem.hpp"
#include "rs_basic.hpp"


using namespace std;
using namespace mfem;

extern int num_equation;

///
/// Implements a simple Rusanov flux
///
class RiemannSolverRusanov : public RiemannSolver
{

public:

    /// Constructor
    RiemannSolverRusanov() : RiemannSolver() {};

    /// Destructor
    ~RiemannSolverRusanov() {};

    /// Compute numerical flux
    virtual double Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, Vector &flux, bool debug = false) override;
};

#endif // RIEMANN_SOLVER_RUS_H
