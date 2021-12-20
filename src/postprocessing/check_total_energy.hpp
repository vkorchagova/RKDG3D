#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

double ComputeTotalEnergy(ParMesh* mesh, ParFiniteElementSpace* vfes, ParGridFunction& sol);

// void ComputeU(ParGridFunction& sol, ParFiniteElementSpace* dfes, ParGridFunction& U);