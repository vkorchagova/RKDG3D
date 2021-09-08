// MFEM Based DG Application fro Solving Discontinuous Gas Dynamics Problems
//
// Compile with: make -j 12
//
// Sample runs:
//
//         dgns -p 1 -r 2 -o 1 -s 3
//
// Case settings could be placed into the YAML file
// Reading by rapidYAML open-source library
//
// Description:  This code solves the compressible Euler system of
//                    equations, a model nonlinear hyperbolic PDE, with a
//                    discontinuous Galerkin (DG) formulation.
//
//                    Note that as the order of the spatial discretization increases,
//                    the timestep must become smaller. This example currently uses a
//                    simple estimate derived by Cockburn and Shu for the 1D RKDG
//                    method. An additional factor can be tuned by passing the --cfl
//                    (or -c shorter) flag.
//                    
//                    Based on MFEM Example 18 for continuous solutions

#include "mfem.hpp"
#include <fstream>
#include <sstream>
#include <iostream>


#include "dg.hpp"
#include "case_manager.hpp"
#include <mpi.h>

// Default value for the problem setup
int problem = 0;

// Equation constant parameters
int num_equation = 4;

// Physics parameters (updated by case)
double specific_heat_ratio = 1.4;
const double gas_constant = 1.0;
double covolume_constant = 0.0;

// Maximum characteristic speed (updated by integrators)
double max_char_speed;

// rank proc
int myRank;

// number of procs
int numProcs;

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
    // x_old.SetFromTrueDofs(u_block_old);
    // x.SetFromTrueDofs(u_block);
    // rhok.SetFromTrueDofs(u_block.GetBlock(0));
    // mom.SetFromTrueDofs(u_block.GetBlock(1));
    // energy.SetFromTrueDofs(u_block.GetBlock(num_equation-1));

    // cout << "IN Update 1 x = " << endl;
    // x.Print(cout);

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
    // x_old.Update();
    x.Update();
    avgr.updateSolutions();


    // cout << "IN Update 2 x = " << endl;
    // x.Print(cout);
   
    // rhok.Update();
    // mom.Update();
    // energy.Update();

    // rhoInd.Update();
    // rhoUInd.Update();
    // rhoVInd.Update();
    
    // if (num_equation == 5) 
    //     rhoWInd.Update();
    // EInd.Update(); 

      


    if (pmesh.Nonconforming())
    {
        // cout << "if (pmesh.Nonconforming())" << endl;
        // Load balance the mesh.
        pmesh.Rebalance();

        // Update the space again, this time a GridFunction redistribution matrix
        // is created. Apply it to the solution.
        fes.Update();
        dfes.Update();
        vfes.Update();
        fes_const.Update();
        avgr.updateSpaces();

        // x_old.Update();
        x.Update();
        
        rhok.Update();
        mom.Update();
        energy.Update();

        // rhoInd.Update();
        // rhoUInd.Update();
        // rhoVInd.Update();
        
        // if (num_equation == 5) 
        //     rhoWInd.Update();
        // EInd.Update();
        avgr.updateSolutions();

        // cout << " end Nonconforming()" << endl;
    } 

    // Compute new offsets
    for (int k = 0; k <= num_equation; k++) 
        offsets[k] = k * vfes.GetNDofs();
    for (int k = 0; k <= num_equation; k++) 
        offsets_const[k] = k * fes_const.GetNDofs();

    u_block.Update(x,offsets);
    u_ind.Update(offsets_const);
    // u_block_old.Update(offsets);

    // cout << "IN Update 3 u_block = " << endl;
    // u_block.Print(cout);

    rhok.MakeRef(&fes,x,offsets[0]);
    mom.MakeRef(&dfes,x,offsets[1]);
    energy.MakeRef(&fes,x,offsets[num_equation-1]);

    // rhok.GetTrueDofs(u_block.GetBlock(0));
    // mom.GetTrueDofs(u_block.GetBlock(1));
    // energy.GetTrueDofs(u_block.GetBlock(num_equation-1));

    // // x_old.GetTrueDofs(u_block_old);
    // x.GetTrueDofs(u_block);

    // cout << "IN Update 4 x = " << endl;
    // x.Print(cout);

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


int main(int argc, char *argv[])
{
    // 0. Initialize MPI.
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    // 1. Prepare for case file parsing

    // 1.1. Set default values of parameters needed in main.cpp 
    bool restart = false;
    int restart_cycle = 0;
    int max_restart_frames = 10;
    double t_final = 1.0;
    double dt = 0.001;
    double cfl = -1;
    bool visualization = false;
    bool paraview = true;
    int vis_steps = 1;
    int precision = 8;
    cout.precision(precision);

    int max_ref_it = 1;

    // 1.2. Read case settings filename
    std::string fileName = "settings.yml";
    
    
    // 2. Read the mesh from the given mesh file.
    VisItDataCollection restart_dc(MPI_COMM_WORLD, "restart");

    CaseManager manager(fileName, restart_dc);

    ParMesh *pmesh;
    manager.loadMesh(pmesh);

    // 2.1. Update num_eqn for 2d or 3d problem

    int dim = pmesh->Dimension();
    int sdim = pmesh->SpaceDimension();
    manager.checkNumEqn(dim);


    // 3. Define the discontinuous DG finite element space of the given
    //     polynomial order on the refined mesh.

    int order = manager.getSpatialOrder();

    DG_FECollection fec(order, dim);
    // Finite element space for a scalar (thermodynamic quantity)
    ParFiniteElementSpace fes(pmesh, &fec);
    // Finite element space for a mesh-dim vector quantity (momentum)
    ParFiniteElementSpace dfes(pmesh, &fec, dim, Ordering::byNODES);
    // Finite element space for all variables together (total thermodynamic state)
    ParFiniteElementSpace vfes(pmesh, &fec, num_equation, Ordering::byNODES);
    // This example depends on this ordering of the space.
    MFEM_ASSERT(fes.GetOrdering() == Ordering::byNODES, "");

    HYPRE_Int glob_size = vfes.GlobalTrueVSize();
    
    if (myRank == 0) 
        cout << "Number of unknowns: " << glob_size << endl; 


    // 4. Define the initial conditions, save the corresponding mesh and grid
    //     functions to a file.

    // The solution u has components {density, x-momentum, y-momentum, (z-momentum), energy}.
    // These are stored contiguously in the BlockVector u_block.
    Array<int> offsets(num_equation + 1);
    for (int k = 0; k <= num_equation; k++) 
        offsets[k] = k * vfes.GetNDofs();
    // BlockVector u_block(offsets);
    // ParGridFunction sol(&vfes, u_block.GetData());
    ParGridFunction sol(&vfes);
    BlockVector u_block(sol, offsets);
    ParGridFunction sol_old(&vfes);
    BlockVector u_block_old(sol_old, offsets);


    // 5. Set up the nonlinear form corresponding to the DG discretization of the
    //     flux divergence, and assemble the corresponding mass matrix.
    MixedBilinearForm Aflux(&dfes, &fes);
    DomainIntegrator* integ = new DomainIntegrator(dim);
    Aflux.AddDomainIntegrator(integ);
    Aflux.Assemble();

    ParNonlinearForm A(&vfes);

    RiemannSolver* rsolver = NULL;
    manager.loadRiemannSolver(rsolver);

    A.AddInteriorFaceIntegrator(new FaceIntegrator(*rsolver, dim));
    manager.addBoundaryIntegrators(A, *pmesh, *rsolver, dim);


    // 6. Define the time-dependent evolution operator describing the ODE
    //     right-hand side, and perform time-integration (looping over the time
    //     iterations, ti, with a time-step dt).
    FE_Evolution euler(vfes, A, &(Aflux.SpMat()));


    // 7. Define limiter and troubled cells indicator

    // 7.1. Define additional block vector for indicator values to have a possibility of visualization
    DG_FECollection fec_const(0, dim);
    ParFiniteElementSpace fes_const(pmesh, &fec_const);
    Array<int> offsets_const(num_equation + 1);
    for (int k = 0; k <= num_equation; k++) { offsets_const[k] = k * fes_const.GetNDofs(); }
    BlockVector indicatorData(offsets_const);

    // 7.2. Load averager, indicator and limiter objects
    Averager avgr(&vfes, offsets, dim);
    Indicator *ind = NULL;
    Limiter *l = NULL;
    manager.loadLimiter(avgr,ind,l,offsets,dim,indicatorData,vfes);

    // 8. Define the ODE solver used for time integration. Several explicit
    //     Runge-Kutta methods are available.


    ODESolver *ode_solver = NULL;
    manager.loadTimeSolver(ode_solver, l);
    manager.loadTimeSettings(t_final, dt, cfl);

    // 9. Determine the minimum element size.

    double hmin = 0.0;
    if (cfl > 0)
    {
        double my_hmin = pmesh->GetElementSize(0, 1);
        for (int i = 1; i < pmesh->GetNE(); i++)
        {
            my_hmin = min(pmesh->GetElementSize(i, 1), my_hmin);
        }
        // Reduce to find the global minimum element size
        MPI_Allreduce(&my_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, pmesh->GetComm());
    }
    
    
    // 10. Create data collection for solution output: either VisItDataCollection for
    // ascii data files, or SidreDataCollection for binary data files.

    manager.getVisSteps(vis_steps);
    // DataCollection *dc = NULL;

    // ParGridFunction rhok(&fes, u_block.GetBlock(0));
    // ParGridFunction mom(&dfes, u_block.GetData() + offsets[1]);
    // ParGridFunction energy(&fes, u_block.GetBlock(dim+1));

    // ParGridFunction rhoInd(&fes_const, indicatorData.GetBlock(0));
    // ParGridFunction rhoUInd(&fes_const, indicatorData.GetBlock(1));
    // ParGridFunction rhoVInd(&fes_const, indicatorData.GetBlock(2));
    // ParGridFunction rhoWInd;
    // if (dim == 3) 
    //     rhoWInd = ParGridFunction(&fes_const, indicatorData.GetBlock(3));
    // ParGridFunction EInd(&fes_const, indicatorData.GetBlock(dim+1));


    ParGridFunction rhok, mom, energy;
    ParGridFunction rhoInd, rhoUInd, rhoVInd, rhoWInd, EInd;

    rhok.MakeRef(&fes, sol, offsets[0]);
    mom.MakeRef(&dfes, sol, offsets[1]);
    energy.MakeRef(&fes, sol, offsets[dim+1]);

    // rhoInd.MakeRef(&fes_const, indicatorData, offsets_const[0]);
    // rhoUInd.MakeRef(&fes_const, indicatorData, offsets_const[1]);
    // rhoVInd.MakeRef(&fes_const, indicatorData, offsets_const[2]);
    // if (dim == 3)
    //     rhoWInd.MakeRef(&fes_const, indicatorData, offsets_const[3]);
    // EInd.MakeRef(&fes_const, indicatorData, offsets_const[dim+1]);


    //////// TRYING TO RUN DYNAMIC MESH REFINEMENT
    DG_FECollection flux_fec(order, dim);
    RT_FECollection smooth_flux_fec(order-1, dim);
    auto flux_fes = &fes;//new ParFiniteElementSpace(pmesh, &flux_fec, sdim);
    auto smooth_flux_fes = new ParFiniteElementSpace(pmesh, &smooth_flux_fec);
    ErrorEstimator* estimator = 0;
    ThresholdRefiner* refiner = 0;
    ThresholdDerefiner* derefiner = 0;

    if (manager.is_adaptive())
    {
        estimator = new L2ZienkiewiczZhuEstimator(*integ, rhok, flux_fes, smooth_flux_fes);
        refiner = new ThresholdRefiner(*estimator);
        derefiner = new ThresholdDerefiner(*estimator);
        manager.loadAdaptiveMeshSettings(*refiner,*derefiner);
    }

    ////////

    // 7.3. LIMIT INITIAL CONDITIONs

    if (myRank == 0) cout << "limit sol OK\n";

    if (manager.is_adaptive())
    {
        refiner->Reset();
        derefiner->Reset();
    }

        // cout << "=== before ref iters === \n" << endl;

    for (int ref_it = 1; ; ref_it++)
    {
            // cout << "--- REF ITER #" << ref_it << endl;
            // cout << "... perform tstep ..." << endl;
            // sol.Print(cout);

        manager.loadInitialSolution(vfes, offsets, u_block, sol);
        pmesh->ExchangeFaceNbrData(); 
        sol.ExchangeFaceNbrData();

        l->update(sol);

            // cout << "... after tstep (elements num = " << pmesh->GetNE() << ")" << endl;
            // sol.Print(cout);

            // cout << "VIS" << endl;
            // rhok.Print(cout);

            // cout << "ublock[0]" << endl;
            // u_block.GetBlock(0).Print(cout);

        if (manager.is_adaptive())
        {
            refiner->Apply(*pmesh);

            // cout << "... after ref (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")" << endl;
            // sol.Print(cout);

            // 22. Update the space, interpolate the solution, rebalance the mesh.
            UpdateAndRebalance(
                *pmesh, 
                fes, 
                dfes, 
                vfes, 
                fes_const, 
                sol, 
                sol_old, 
                Aflux, 
                A, 
                rhok, 
                mom, 
                energy, 
                rhoInd, 
                rhoUInd, 
                rhoVInd, 
                rhoWInd, 
                EInd,
                u_block,
                u_block_old,
                indicatorData,
                offsets,
                offsets_const,
                avgr
            );

            // cout << "... after rebalance " << vfes.GlobalTrueVSize() << endl;
            // sol.Print(cout);

            euler.UpdateAfluxPointer(&(Aflux.SpMat()));
            euler.UpdateInverseMassMatrix();

            // sol.Print(cout);

            if (refiner->Stop())
            {
                break;
            }
            // cout << "... after sol = sol_old " << endl; 
        }
        else
        {
            break;
        }
    }
        // cout << "... before deref" << endl;
    if (manager.is_adaptive())
    {
        if (derefiner->Apply(*pmesh))
        {
            if (myRank == 0)
            {
                // cout << "\nDerefined elements." << endl;
                cout << "... after deref (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")" << endl;
            }

            // 24. Update the space and the solution, rebalance the mesh.
            UpdateAndRebalance(
                *pmesh, 
                fes, 
                dfes, 
                vfes, 
                fes_const, 
                sol, 
                sol_old, 
                Aflux, 
                A, 
                rhok, 
                mom, 
                energy, 
                rhoInd, 
                rhoUInd, 
                rhoVInd, 
                rhoWInd, 
                EInd,
                u_block,
                u_block_old,
                indicatorData,
                offsets,
                offsets_const,
                avgr
            );

            // cout << "after second basic rebalance" << endl;

            euler.UpdateAfluxPointer(&(Aflux.SpMat()));
            euler.UpdateInverseMassMatrix();

            // cout << "after second rebalance" << endl;
            // sol.Print(cout);
        }

        // UpdateAndRebalance(
        //         *pmesh, 
        //         fes, 
        //         dfes, 
        //         vfes, 
        //         fes_const, 
        //         sol, 
        //         sol_old, 
        //         Aflux, 
        //         A, 
        //         rhok, 
        //         mom, 
        //         energy, 
        //         rhoInd, 
        //         rhoUInd, 
        //         rhoVInd, 
        //         rhoWInd, 
        //         EInd,
        //         u_block,
        //         u_block_old,
        //         indicatorData,
        //         offsets,
        //         offsets_const,
        //         avgr
        //     );
        // // cout << "after last basic rebalance" << endl;           
        // euler.UpdateAfluxPointer(&(Aflux.SpMat()));
        // euler.UpdateInverseMassMatrix(); 
        cout << "after last rebalance (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")"  << endl;
    } // end adaptive mesh

    l->update(sol);



    ParaViewDataCollection *pd = NULL;
    if (paraview)
    {
        pd = new ParaViewDataCollection("PV", pmesh);
        pd->RegisterField("mom", &mom);
        pd->RegisterField("rho", &rhok);
        pd->RegisterField("energy", &energy);
        // pd->RegisterField("rhoInd", &rhoInd);
        // pd->RegisterField("rhoUInd", &rhoUInd);
        // pd->RegisterField("rhoVInd", &rhoVInd);
        // if (dim == 3) pd->RegisterField("rhoWInd", &rhoWInd);
        // pd->RegisterField("EInd", &EInd);

        pd->SetLevelsOfDetail(1);
        pd->SetCycle(0);
        pd->SetTime(0.0);
        pd->Save();
    }


    // 11. Initialize restart queue for control of number of saved frames
    manager.initializeRestartQueue();

    
    // 12. Start the timer.
    
    tic_toc.Clear();
    tic_toc.Start();

    double t = manager.setStartTime();
    euler.SetTime(t);
    ode_solver->Init(euler);

    socketstream sout;

    if (!manager.is_restart())
    {
        if (cfl > 0)
        {
            // Find a safe dt, using a temporary vector. Calling Mult() computes the
            // maximum char speed at all quadrature points on all faces.
            max_char_speed = 0.;
            Vector z(sol.Size());
            A.Mult(sol, z);
            // Reduce to find the global maximum wave speed
            {
                double all_max_char_speed;
                MPI_Allreduce(&max_char_speed, &all_max_char_speed,
                                  1, MPI_DOUBLE, MPI_MAX, pmesh->GetComm());
                max_char_speed = all_max_char_speed;
            }
            dt = cfl * hmin / max_char_speed / (2*order+1);
        }
    }
    else
    {
        dt = restart_dc.GetTimeStep();
    }


    // 13. Integrate in time.

    if (myRank == 0) cout << "START TIME CYCLE" << endl;

    bool done = false;
    double t_real = 0.0;

    // sol_old = sol;

    for (int ti = manager.setStartTimeCycle(); !done; )
    {
        t_real = tic_toc.RealTime();
        double dt_real = min(dt, t_final - t);

        // Make sure errors will be recomputed in the following.
        if (manager.is_adaptive())
        {
            refiner->Reset();
            derefiner->Reset();
        }

        // cout << "=== before ref iters === \n" << endl;

        for (int ref_it = 1; ref_it <= max_ref_it; ref_it++)
        {
            // cout << "--- REF ITER #" << ref_it << endl;
            // cout << "... perform tstep ..." << endl;
            // sol.Print(cout);

            ode_solver->Step(sol, t, dt_real);

            // cout << "... after tstep (elements num = " << pmesh->GetNE() << ")" << endl;
            // sol.Print(cout);

            // cout << "VIS" << endl;
            // rhok.Print(cout);

            // cout << "ublock[0]" << endl;
            // u_block.GetBlock(0).Print(cout);

            if (manager.is_adaptive())
            {
                refiner->Apply(*pmesh);

                // cout << "... after ref (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")" << endl;
                // sol.Print(cout);

                // 22. Update the space, interpolate the solution, rebalance the mesh.
                UpdateAndRebalance(
                    *pmesh, 
                    fes, 
                    dfes, 
                    vfes, 
                    fes_const, 
                    sol, 
                    sol_old, 
                    Aflux, 
                    A, 
                    rhok, 
                    mom, 
                    energy, 
                    rhoInd, 
                    rhoUInd, 
                    rhoVInd, 
                    rhoWInd, 
                    EInd,
                    u_block,
                    u_block_old,
                    indicatorData,
                    offsets,
                    offsets_const,
                    avgr
                );

                // cout << "... after rebalance " << vfes.GlobalTrueVSize() << endl;
                // sol.Print(cout);

                euler.UpdateAfluxPointer(&(Aflux.SpMat()));
                euler.UpdateInverseMassMatrix();

                // sol.Print(cout);

                if (refiner->Stop() || ref_it == max_ref_it)
                {
                    // Aflux.Update(); // Free the assembled data
                    // A.Update();
                    // cout << "... stop refiner" << endl;
                    break;
                }

                // load to sol old solution again just to make the same time step
                sol = sol_old;
                // cout << "... after sol = sol_old " << endl; 
            }
            else
            {
                break;
            }
        }
        // cout << "... before deref" << endl;
        if (manager.is_adaptive())
        {
            if (derefiner->Apply(*pmesh))
            {
                if (myRank == 0)
                {
                    // cout << "\nDerefined elements." << endl;
                    cout << "... after deref (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")" << endl;
                }

                // 24. Update the space and the solution, rebalance the mesh.
                UpdateAndRebalance(
                    *pmesh, 
                    fes, 
                    dfes, 
                    vfes, 
                    fes_const, 
                    sol, 
                    sol_old, 
                    Aflux, 
                    A, 
                    rhok, 
                    mom, 
                    energy, 
                    rhoInd, 
                    rhoUInd, 
                    rhoVInd, 
                    rhoWInd, 
                    EInd,
                    u_block,
                    u_block_old,
                    indicatorData,
                    offsets,
                    offsets_const,
                    avgr
                );

                // cout << "after second basic rebalance" << endl;

                euler.UpdateAfluxPointer(&(Aflux.SpMat()));
                euler.UpdateInverseMassMatrix();

                // cout << "after second rebalance" << endl;
                // sol.Print(cout);
            }

            // UpdateAndRebalance(
            //         *pmesh, 
            //         fes, 
            //         dfes, 
            //         vfes, 
            //         fes_const, 
            //         sol, 
            //         sol_old, 
            //         Aflux, 
            //         A, 
            //         rhok, 
            //         mom, 
            //         energy, 
            //         rhoInd, 
            //         rhoUInd, 
            //         rhoVInd, 
            //         rhoWInd, 
            //         EInd,
            //         u_block,
            //         u_block_old,
            //         indicatorData,
            //         offsets,
            //         offsets_const,
            //         avgr
            //     );
            // // cout << "after last basic rebalance" << endl;           
            // euler.UpdateAfluxPointer(&(Aflux.SpMat()));
            // euler.UpdateInverseMassMatrix(); 
            cout << "after last rebalance (elements num = "<< pmesh->GetNE() << ";" << vfes.GlobalTrueVSize() << ")"  << endl;
            //     sol.Print(cout);


            // sol_old just for making steps in refinement
            sol_old = sol;
        } // end adaptive mesh

        if (cfl > 0)
        {
            // Reduce to find the global maximum wave speed
            {
                double all_max_char_speed;
                MPI_Allreduce(&max_char_speed, &all_max_char_speed,
                                  1, MPI_DOUBLE, MPI_MAX, pmesh->GetComm());
                max_char_speed = all_max_char_speed;
            }
            dt = cfl * hmin / max_char_speed / (2*order+1);
        }
        

        ti++;

        t_real = tic_toc.RealTime() - t_real;

        done = (t >= t_final - 1e-8*dt);

        if (done || ti % vis_steps == 0)
        {
            
            MPI_Barrier(pmesh->GetComm());

            if (paraview)
            {
                // rhok.SetFromTrueDofs(u_block.GetBlock(0));
                // mom.SetFromTrueDofs(u_block.GetBlock(1));
                // energy.SetFromTrueDofs(u_block.GetBlock(2));
          
                // rhoInd.SetFromTrueDofs(u_block.GetBlock(0));
                // rhoUInd.SetFromTrueDofs(u_block.GetBlock(1));
                // rhoVInd.SetFromTrueDofs(u_block.GetBlock(2));
                // if (dim == 3)
                //     rhoWInd.SetFromTrueDofs(u_block.GetBlock(3));
                // EInd.SetFromTrueDofs(u_block.GetBlock(dim+1));

                //  rhoInd.MakeTRef(&fes_const, indicatorData, offsets_const[0]);
                // rhoUInd.MakeTRef(&fes_const, indicatorData, offsets_const[1]);
                // rhoVInd.MakeTRef(&fes_const, indicatorData, offsets_const[2]);
                // if (dim == 3)
                //     rhoWInd.MakeTRef(&fes_const, indicatorData, offsets_const[3]);
                // EInd.MakeTRef(&fes_
                
                pd->SetCycle(ti);
                pd->SetTime(t);
                pd->Save();
                if (myRank == 0) {cout << "ParaView OK\n";}
                
            }

            restart_dc.SetCycle(ti);
            restart_dc.SetTime(t);
            restart_dc.SetTimeStep(dt);
            restart_dc.Save();

            if (myRank == 0) 
            {
                manager.cleanPreviousRestartFrames(ti);
            }
            
        }
        

        if (myRank == 0)
        {
            cout << "Time step: " << ti << "\tdt: " << dt << "\tPhys time: " << t << " \tExecution time: " << tic_toc.RealTime() << " \tReal dt time: " << t_real << endl;
        }
    }

    tic_toc.Stop();

    if (myRank == 0) 
    { 
        cout << " done, " << tic_toc.RealTime() << "s." << endl; 
    }

    delete estimator;

    // Free the used memory.
    delete derefiner;
    delete refiner;
    delete estimator;
    delete smooth_flux_fes;

    delete pd;
    delete ode_solver;
    delete l;
    delete ind;
    // delete rsolver;
    // delete integ;
    // delete pmesh;
    cout << "before mpi finalize" << endl;
    MPI_Finalize();
    cout << "after mpi finalize" << endl;

    return 0;
}


