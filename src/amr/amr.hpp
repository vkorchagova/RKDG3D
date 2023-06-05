#include "mfem.hpp"
#include "averager.hpp"

using namespace mfem;

/// Number of equations
extern int num_equation;


// class AMRController
// {

//    /// Pointer to mesh
//    ParMesh* pmesh;

//    /// List of spaces to be updated
//    List<ParFiniteElementSpace&> spaces;

   

// public:

//    void UpdateAndRebalance(
//       ParMesh &pmesh,
//       List<ParFiniteElementSpace&> spaces,
//       List<ParGridFunction&> fields,
//       MixedBilinearForm &a, /* Aflux */
//       ParNonlinearForm &b,  /* A */
//       BlockVector &u_block,
//       BlockVector &u_block_old,
//       Array<int> &offsets,
//       Array<int> &offsets_const,
//       Averager& avgr
//    )
// }



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
   BlockVector &u_block,
   BlockVector &u_block_old,
   ParGridFunction &u_ind,
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