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


#include "fe_evolution.hpp"
#include "domain_integrator.hpp"
#include "face_integrator.hpp"
#include "case_manager.hpp"
#include "amr.hpp"
#include "check_total_energy.hpp"
#include <mpi.h>

// Default value for the problem setup
int problem = 0;

// Equation constant parameters
int num_equation = 4;

// Physics parameters (updated by case)
double specific_heat_ratio = 1.4;
double gas_constant = 287;
double covolume_constant = 0.0;

// Maximum characteristic speed (updated by integrators)
double max_char_speed;

// rank proc
int myRank;

// number of procs
int numProcs;


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

    if (myRank == 0) 
    {
        cout << "Dimension: " << dim << endl;
        cout << "Number of equations: " << num_equation << endl;
    }
    

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
    ParGridFunction indicatorData(&fes_const);

    // 7.2. Load averager, indicator and limiter objects
    Averager avgr(&fes, offsets, dim);
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
    
    ParGridFunction rhok, mom, energy;
    ParGridFunction rhoInd;

    rhok.MakeRef(&fes, sol, offsets[0]);
    mom.MakeRef(&dfes, sol, offsets[1]);
    energy.MakeRef(&fes, sol, offsets[dim+1]);

    ParGridFunction U(&dfes);
    ParGridFunction p(&fes);
    ParGridFunction T(&fes);

    ParGridFunction UMean(&dfes);
    ParGridFunction pMean(&fes);
    ParGridFunction TMean(&fes);
    ParGridFunction rhoMean(&fes);

    // 11. Initialize spaces and objects for dynamic mesh refinement

    DG_FECollection flux_fec(order, dim);
    RT_FECollection smooth_flux_fec(max(order-1,0), dim);
    ParFiniteElementSpace* flux_fes = &dfes;
    ParFiniteElementSpace* smooth_flux_fes = new ParFiniteElementSpace(pmesh, &smooth_flux_fec);
    ErrorEstimator* estimator = 0;
    ThresholdRefiner* refiner = 0;
    ThresholdDerefiner* derefiner = 0;

    if (manager.is_adaptive())
    {
        estimator = new L2ZienkiewiczZhuEstimator(*integ, rhok, *flux_fes, *smooth_flux_fes);
        refiner = new ThresholdRefiner(*estimator);
        derefiner = new ThresholdDerefiner(*estimator);
        manager.loadAdaptiveMeshSettings(*refiner,*derefiner);
    }

    // 12. Limit initial conditions and refine mesh accordingly to them

    if (manager.is_adaptive())
    {
        refiner->Reset();
        derefiner->Reset();
    }


    for (int ref_it = 1; ; ref_it++)
    {
        manager.loadInitialSolution(vfes, offsets, u_block, sol);
        l->update(sol);

        if (myRank == 0) cout << "Refinement iteration # " << ref_it << "..." << endl;

        if (manager.is_adaptive() && !manager.is_restart())
        {
            refiner->Apply(*pmesh);

            // Update the space, interpolate the solution, rebalance the mesh.
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
                u_block,
                u_block_old,
                indicatorData,
                offsets,
                offsets_const,
                avgr,
                U,
                p,
                T,
                UMean,
                pMean,
                TMean,
                rhoMean
            );

            // cout << "... after rebalance " << vfes.GlobalTrueVSize() << endl;
            // sol.Print(cout);

            euler.UpdateAfluxPointer(&(Aflux.SpMat()));
            euler.UpdateInverseMassMatrix();

            if (refiner->Stop() )
            {
                break;
            }
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
                cout << "... after deref (elements num = "<< pmesh->GetNE() << ";" << vfes.TrueVSize() << ";"
                     << vfes.GlobalTrueVSize() << ")" << endl;
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
                u_block,
                u_block_old,
                indicatorData,
                offsets,
                offsets_const,
                avgr,
                U,
                p,
                T,
                UMean,
                pMean,
                TMean,
                rhoMean
            );

            // cout << "after second basic rebalance" << endl;

            euler.UpdateAfluxPointer(&(Aflux.SpMat()));
            euler.UpdateInverseMassMatrix();

            // cout << "after second rebalance" << endl;
            // sol.Print(cout);
        }

        l->update(sol);

        cout << "after last rebalance (elements num = "
             << pmesh->GetNE() << ";" 
             << vfes.TrueVSize() << ";"
             << vfes.GlobalTrueVSize() << ")"  << endl;

        sol_old = sol;
    } // end adaptive mesh

    Vector localState(num_equation);

    int fesNDofs = fes.GetNDofs();

    for (int i = 0; i < fesNDofs; ++i)
    {
        U[i] = mom[i]/rhok[i];
        U[i+fesNDofs] = mom[i+fesNDofs]/rhok[i];
        if (dim == 3) U[i+2*fesNDofs] = mom[i+2*fesNDofs]/rhok[i];

        localState[0] = rhok[i];
        localState[1] = mom[i];
        localState[2] = mom[i+fesNDofs];
        if (dim == 3) localState[3] = mom[i+2*fesNDofs];
        localState[num_equation-1] = energy[i];

        p[i] = ComputePressure(localState, dim);
        T[i] = ComputeTemperature(localState, dim);

        UMean[i] = 0.0;
        UMean[i+fesNDofs] = 0.0;
        if (dim == 3) UMean[i+2*fesNDofs] = 0.0;
        pMean[i] = 0.0;
        TMean[i] = 0.0;
        rhoMean[i] = 0.0;
    }


    if (myRank == 0) cout << "Projection of initial conditions OK" << endl;

    // 13. Initialize data collection for ParaView
    ParaViewDataCollection *pd = NULL;
    if (paraview)
    {
        pd = new ParaViewDataCollection("PV", pmesh);
        pd->RegisterField("mom", &mom);
        pd->RegisterField("rho", &rhok);
        pd->RegisterField("energy", &energy);

        if (manager.write_indicators())
        {
            pd->RegisterField("indicator", &indicatorData);
        }

        pd->RegisterField("U", &U);
        pd->RegisterField("p", &p);
        pd->RegisterField("T", &T);

        pd->RegisterField("UMean", &UMean);
        pd->RegisterField("pMean", &pMean);
        pd->RegisterField("TMean", &TMean);
        pd->RegisterField("rhoMean", &rhoMean);

        pd->SetLevelsOfDetail(1);
        pd->SetCycle(0);
        pd->SetTime(0.0);
        pd->Save();

        if (myRank == 0) cout << "Paraview OK" << endl;
    }

    // 14. Initialize restart queue for control of number of saved frames
    manager.initializeRestartQueue();

    // 15. Start the timer.
    
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


    restart_dc.SetCycle(0);
    restart_dc.SetTime(0);
    restart_dc.SetTimeStep(dt);
    restart_dc.Save();

    // 16. Integrate in time.

    if (myRank == 0) cout << "START TIME CYCLE" << endl;

    bool done = false;
    double t_real = 0.0;
    double e_tot = 0.0;


    Array<int> vdofs;
    Vector el_x;
    
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


        for (int ref_it = 1; ; ref_it++)
        {
            ode_solver->Step(sol, t, dt_real);

            // vfes.GetElementVDofs(50, vdofs);
            // sol.GetSubVector(vdofs, el_x);

            // el_x.Print(cout << std::setprecision(30) << "sol 50 = ");

            // vfes.GetElementVDofs(1562, vdofs);
            // sol.GetSubVector(vdofs, el_x);

            // el_x.Print(cout << std::setprecision(30) << "sol 1562 = ");

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

                // if (myRank == 0)
                // {
                //     cout << "hmin = " << hmin << "; max_char_speed = " << max_char_speed << "; 2*order+1 = " << 2*order+1 << endl;
                // }
            }

            // cout << "--- REF ITER #" << ref_it << endl;
            // cout << "... perform tstep for rank = " << myRank << "; (elements num = " << pmesh->GetNE() << ")" << endl;

            // cout << "VIS" << endl;
            // rhok.Print(cout);

            // cout << "ublock[0]" << endl;
            // u_block.GetBlock(0).Print(cout);

            if (manager.is_adaptive())
            {
                if (myRank == 0) cout << "Refinement iteration # " << ref_it << "..." << endl;
            
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
                    u_block,
                    u_block_old,
                    indicatorData,
                    offsets,
                    offsets_const,
                    avgr,
                    U,
                    p,
                    T,
                    UMean,
                    pMean,
                    TMean,
                    rhoMean
                );

                // cout << "... after rebalance " << vfes.GlobalTrueVSize() << endl;
                // sol.Print(cout);

                euler.UpdateAfluxPointer(&(Aflux.SpMat()));
                euler.UpdateInverseMassMatrix();

                // sol.Print(cout);

                if (refiner->Stop())
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
                    u_block,
                    u_block_old,
                    indicatorData,
                    offsets,
                    offsets_const,
                    avgr,
                    U,
                    p,
                    T,
                    UMean,
                    pMean,
                    TMean,
                    rhoMean
                );

                // cout << "after second basic rebalance" << endl;

                euler.UpdateAfluxPointer(&(Aflux.SpMat()));
                euler.UpdateInverseMassMatrix();

                // cout << "after second rebalance" << endl;
                // sol.Print(cout);
            }

            
            cout << "after last rebalance (elements num = "
             << pmesh->GetNE() << ";" 
             << vfes.TrueVSize() << ";"
             << vfes.GlobalTrueVSize() << ")"  << endl;

            l->update(sol);

            // sol_old just for making steps in refinement
            sol_old = sol;
        } // end adaptive mesh



        // update time step
        t += dt;

        ti++;

        if (manager.check_total_energy())
            e_tot = ComputeTotalEnergy(pmesh, &vfes, sol);

        t_real = tic_toc.RealTime() - t_real;

        for (int i = 0; i < fes.GetNDofs(); ++i)
        {
            U[i] = mom[i]/rhok[i];
            U[i+fes.GetNDofs()] = mom[i+fes.GetNDofs()]/rhok[i];
            if (dim == 3) U[i+2*fes.GetNDofs()] = mom[i+2*fes.GetNDofs()]/rhok[i];

            localState[0] = rhok[i];
            localState[1] = mom[i];
            localState[2] = mom[i+fes.GetNDofs()];
            if (dim == 3) localState[3] = mom[i+2*fes.GetNDofs()];
            localState[num_equation-1] = energy[i];

            p[i] = ComputePressure(localState, dim);
            T[i] = ComputeTemperature(localState, dim);

            UMean[i] += U[i] * dt;
            UMean[i+fesNDofs] += U[i+fesNDofs] * dt;
            if (dim == 3) UMean[i+2*fesNDofs] += U[i+2*fesNDofs] * dt;
            pMean[i] += p[i] * dt;
            TMean[i] += T[i] * dt;
            rhoMean[i] += rhok[i] * dt;
        }

        done = (t >= t_final - 1e-8*dt);

        if (done || ti % vis_steps == 0)
        {   
            MPI_Barrier(pmesh->GetComm());

            for (int i = 0; i < fes.GetNDofs(); ++i)
            {
                UMean[i] /= t;
                UMean[i+fesNDofs] /= t;
                if (dim == 3) UMean[i+2*fesNDofs] /= t;
                pMean[i] /= t;
                TMean[i] /= t;
                rhoMean[i] /= t;
            }

            if (paraview)
            {
                pd->SetCycle(ti);
                pd->SetTime(t);
                pd->Save();
                if (myRank == 0) {cout << "ParaView OK\n";}
            }

            restart_dc.SetCycle(ti);
            restart_dc.SetTime(t);
            restart_dc.SetTimeStep(dt);
            restart_dc.Save();

            for (int i = 0; i < fes.GetNDofs(); ++i)
            {
                UMean[i] *= t;
                UMean[i+fesNDofs] *= t;
                if (dim == 3) UMean[i+2*fesNDofs] *= t;
                pMean[i] *= t;
                TMean[i] *= t;
                rhoMean[i] *= t;
            }

            if (myRank == 0) 
            {
                manager.cleanPreviousRestartFrames(ti);
            }
            
        }
        
        if (myRank == 0)
        {
            cout << "Time step: " << ti 
                 << "\tdt: " << dt 
                 << "\tPhys time: " << t 
                 << " \tExecution time: " << tic_toc.RealTime() 
                 << " \tReal dt time: " << t_real 
                 << " \t E_total: " << setprecision(18) << e_tot /*(manager.check_total_energy() ? e_tot : "NOT COMPUTED");*/
                 << setprecision(6)
                 << endl;
        }

    }

    tic_toc.Stop();

    if (myRank == 0) 
    { 
        cout << " done, " << tic_toc.RealTime() << "s." << endl; 
    }

    // Free the used memory.
    delete pd;
    delete derefiner;
    delete refiner;
    delete estimator; 
    delete smooth_flux_fes;   
    delete ode_solver;
    delete l;
    delete ind;
    delete rsolver;
    delete pmesh;
    MPI_Finalize();

    return 0;
}


