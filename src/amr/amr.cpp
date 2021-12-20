#include "amr.hpp"

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
    Averager& avgr
)
{
    // Update the space: recalculate the number of DOFs and construct a matrix
    // that will adjust any GridFunctions to the new mesh state.
    
    fes.Update(); //false in ()
    dfes.Update();
    vfes.Update();
    fes_const.Update();
    avgr.updateSpaces();

    // Interpolate the solution on the new mesh by applying the transformation
    // matrix computed in the finite element space. Multiple GridFunctions could
    // be updated here.
    x.Update();
    x_old.Update();

    // Compute new offsets
    for (int k = 0; k <= num_equation; k++) 
        offsets[k] = k * vfes.GetNDofs();
    for (int k = 0; k <= num_equation; k++) 
        offsets_const[k] = k * fes_const.GetNDofs();

    u_block.Update(x,offsets);
    u_ind.Update(offsets_const);

    rhok.MakeRef(&fes,x,offsets[0]);
    mom.MakeRef(&dfes,x,offsets[1]);
    energy.MakeRef(&fes,x,offsets[num_equation-1]);
    avgr.updateSolutions();

    if (pmesh.Nonconforming())
    {
        // cout << "if (pmesh.Nonconforming())" << endl;
        // Load balance the mesh.
        pmesh.Rebalance();
        pmesh.ExchangeFaceNbrData(); 

        // Update the space again, this time a GridFunction redistribution matrix
        // is created. Apply it to the solution.
        fes.Update();
        dfes.Update();
        vfes.Update();
        fes_const.Update();
        avgr.updateSpaces();

        x.Update();
        x_old.Update();
        x.ExchangeFaceNbrData();
        x_old.ExchangeFaceNbrData();

        // Compute new offsets
        for (int k = 0; k <= num_equation; k++) 
            offsets[k] = k * vfes.GetNDofs();
        for (int k = 0; k <= num_equation; k++) 
            offsets_const[k] = k * fes_const.GetNDofs();

        u_block.Update(x,offsets);
        u_ind.Update(offsets_const);

        rhok.MakeRef(&fes,x,offsets[0]);
        mom.MakeRef(&dfes,x,offsets[1]);
        energy.MakeRef(&fes,x,offsets[num_equation-1]);

        avgr.updateSolutions();

        // cout << " end Nonconforming()" << endl;
    } 

    // Inform the nonlinear and bilinear forms that the space has changed.
    a.Update();
    b.Update();
    a.Assemble();

    // Free any transformation matrices to save memory.
    fes.UpdatesFinished();
    dfes.UpdatesFinished();
    vfes.UpdatesFinished();
    fes_const.UpdatesFinished();
    avgr.updateFinished();

    // cout << "UpdateAndRebalance OK" << endl;
}