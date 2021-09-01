#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

///
/// Abstract class for Riemann solvers
///
class RiemannSolver
{
public:
    
    /// Physical flux F, left side
    Vector flux1;

    /// Physical flux F, right side
    Vector flux2;

    void Rotate(Vector& state, const Vector& nor, int dim);

    void InverseRotate(Vector& state, const Vector& nor, int dim);

public:

    /// Constructor
    RiemannSolver();

    /// Destructor
    virtual ~RiemannSolver() {};

    /// Compute numerical flux
    virtual double Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, Vector &flux, bool debug = false) = 0;
};

#endif // RIEMANN_SOLVER_H
