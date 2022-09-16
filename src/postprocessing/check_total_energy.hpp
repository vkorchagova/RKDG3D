#include "mfem.hpp"

using namespace mfem;

/// Number of equations
extern int num_equation;

/// Check global conservativity by computation of total energy in system
double ComputeTotalEnergy(ParMesh* mesh, ParFiniteElementSpace* vfes, ParGridFunction& sol);

// void ComputeU(ParGridFunction& sol, ParFiniteElementSpace* dfes, ParGridFunction& U);