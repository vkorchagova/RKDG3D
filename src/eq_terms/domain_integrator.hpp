#ifndef DOMAIN_INTEGRATOR_H
#define DOMAIN_INTEGRATOR_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern const int num_equation;

// Constant (in time) mixed bilinear form multiplying the flux grid function.
// The form is (vec(v), grad(w)) where the trial space = vector L2 space (mesh
// dim) and test space = scalar L2 space.
class DomainIntegrator : public BilinearFormIntegrator
{
private:

    /// Values of shape function at the given point
    Vector shape;

    /// NOT USED
    DenseMatrix flux;

    /// Gradient of test function
    DenseMatrix dshapedr;

    /// Gradient of test function ?
    DenseMatrix dshapedx;

public:
    /// Constructor
    DomainIntegrator(const int dim);

    /// Compute element term with gradients
    virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Tr,
                                       DenseMatrix &elmat);
};

#endif // DOMAIN_INTEGRATOR_H