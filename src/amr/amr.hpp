#include "mfem.hpp"
#include "averager.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

/// Update and rebalance all FEspaces and grid functions after each change of mesh 
void UpdateAndRebalance(
    ParMesh &pmesh, 
    ParFiniteElementSpace &fes,
    ParFiniteElementSpace &dfes, 
    ParFiniteElementSpace &vfes, 
    ParFiniteElementSpace &fes_const,
    ParGridFunction &x, 
    ParGridFunction &x_old, 
    MixedBilinearForm &a, /* Aflux */
    ParNonlinearForm &b,  /* A */
    ParGridFunction &rhok,
    ParGridFunction &mom,
    ParGridFunction &energy,
    ParGridFunction &rhoInd,
    ParGridFunction &rhoUInd,
    ParGridFunction &rhoVInd,
    ParGridFunction &rhoWInd,
    ParGridFunction &EInd,
    BlockVector &u_block,
    BlockVector &u_block_old,
    BlockVector &u_ind,
    Array<int> &offsets,
    Array<int> &offsets_const,
    Averager& avgr,
    ParGridFunction &U,
    ParGridFunction &p,
    ParGridFunction &T,
    ParGridFunction &UMean,
    ParGridFunction &pMean,
    ParGridFunction &TMean,
    ParGridFunction &rhoMean
);