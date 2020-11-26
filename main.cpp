//                                MFEM Example 18 Based DGNS Application
//
// Compile with: make -j12
//
// Sample runs:
//
//       dgns -p 1 -r 2 -o 1 -s 3
//
// Description:  This example code solves the compressible Euler system of
//               equations, a model nonlinear hyperbolic PDE, with a
//               discontinuous Galerkin (DG) formulation.
//
//               Since all
//               boundaries are periodic here, the method's accuracy can be
//               assessed by measuring the difference between the solution and
//               the initial condition at a later time when the vortex returns
//               to its initial location.
//
//               Note that as the order of the spatial discretization increases,
//               the timestep must become smaller. This example currently uses a
//               simple estimate derived by Cockburn and Shu for the 1D RKDG
//               method. An additional factor can be tuned by passing the --cfl
//               (or -c shorter) flag.
//
//               The example demonstrates user-defined bilinear and nonlinear
//               form integrators for systems of equations that are defined with
//               block vectors, and how these are used with an operator for
//               explicit time integrators. In this case the system also
//               involves an external approximate Riemann solver for the DG
//               interface flux. It also demonstrates how to use GLVis and Paraview for
//               in-situ visualization of vector grid functions.
//
//               We recommend viewing examples 9, 14 and 17 before viewing this
//               example.

#include "mfem.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

// Classes FE_Evolution, RiemannSolver, DomainIntegrator and FaceIntegrator
// shared between the serial and parallel version of the example.
#include "dg.hpp"
#include <mpi.h>

// Choice for the problem setup. See InitialCondition in ex18.hpp.
int problem;

// Equation constant parameters
const int num_equation = 4;
double specific_heat_ratio = 1.4;
const double gas_constant = 1.0;
double covolume_constant = 0.0;

// Maximum characteristic speed (updated by integrators)
double max_char_speed;

// rank proc
int myRank;

int numProcs;

int main(int argc, char *argv[])
{
   // 0. Initialize MPI.
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // 1. Parse command-line options.
   problem = 1;
   bool restart = false;
   int restart_cycle = 0;
   const char *mesh_file = "sod2d.msh";
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 1;
   int ode_solver_type = 2;
   double t_final = 0.25;
   double dt = 0.001;
   double cfl = -1;
   bool visualization = false;
   bool visit = false;
   bool paraview = false;
   int vis_steps = 250;
   int indicator_type = 3;
   int limiter_type = 2;
   int riemann_solver_type = 3;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See options in velocity_function().");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly before parallel"
                  " partitioning, -1 for auto.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly after parallel"
                  " partitioning.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Forward Euler,\n\t"
                  "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step. Positive number skips CFL timestep calculation.");
   args.AddOption(&cfl, "-c", "--cfl-number",
                  "CFL number for timestep calculation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView (paraview.org) visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&indicator_type, "-it", "--indicator-type",
                  "Set indicator type: 0 - Nowhere, 1 - Everywhere, 2 - Barth-Jespersen, 3 - Shu");
   args.AddOption(&limiter_type, "-lt", "--limiter-type",
                  "Set limiter type: 1 - FinDiff, 2 - Multiplier.");
   args.AddOption(&riemann_solver_type, "-rst", "--riemann-solver-type",
                  "Set Riemann solver type: 0 - Rusanov, 1 - LLF, 2 - HLL, 3 - HLLC.");
   args.AddOption(&restart, "-restart", "--restart", "-no-restart", "--no-restart",
                  "Enable or disable restart.");
   args.AddOption(&restart_cycle, "-r_cycle", "--restart_cycle",
                  "Value (int) of time cycle for restart.");
   args.Parse();
   if (!args.Good())
   {
      if (myRank == 0) { args.PrintUsage(cout); }
      return 1;
   }
   if (myRank == 0) { args.PrintOptions(cout); }

   // 2. Read the mesh from the given mesh file. This example requires a 2D
   //    periodic mesh, such as ../data/periodic-square.mesh.


   // update num_eqn for 2d or 3d problem
   // if (dim == 2)
   // {
   //    int* ptr = const_cast<int*>(&num_equation);
   //    *ptr = 4;
   // }
   // else if (dim == 3)
   // {
   //    int* ptr =  const_cast<int*>(&num_equation);
   //    *ptr = 5;
   // }

   
   ParMesh *pmesh;
   VisItDataCollection restart_dc(MPI_COMM_WORLD, "restart");

   if (restart) // load par mesh in case of restart
   {
      // 2.1. Read the given data collection.
      restart_dc.Load(restart_cycle);
      if (restart_dc.Error())
      {
         if (myRank == 0)
         {
            cout << "Error loading data collection" << endl;
         }
         return 1;
      }
      // 2.2. Get par mesh pieces from the loaded collection
      pmesh = dynamic_cast<ParMesh*>(restart_dc.GetMesh());
      if (!pmesh)
      {
         if (myRank == 0)
         {
            cout << "The given data collection does not have a parallel mesh."
                 << endl;
         }
         return 2;
      }
   }
   else //load mesh from file and make a partition
   {
      // 2.1. Read the serial mesh on all processors, refine it in serial, then
      //    partition it across all processors and refine it in parallel.
      Mesh *mesh = new Mesh(mesh_file, 1, 1);

      for (int l = 0; l < ser_ref_levels; l++)
      {
         mesh->UniformRefinement();
      }

      if (myRank == 0) { cout << "Number of cells: " << mesh->GetNE() << endl; }

      // 2.2. Define a parallel mesh by a partitioning of the serial mesh. Refine
      //    this mesh further in parallel to increase the resolution. Once the
      //    parallel mesh is defined, the serial mesh can be deleted.
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
      delete mesh;

      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
      restart_dc.SetFormat(1);
      restart_dc.SetMesh(pmesh);  
   }  

   const int dim = pmesh->Dimension();


   // 3. Define the discontinuous DG finite element space of the given
   //    polynomial order on the refined mesh.
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
   if (myRank == 0) { cout << "Number of unknowns: " << glob_size << endl; }

   // 4. Define the initial conditions, save the corresponding mesh and grid
   //    functions to a file.

   // The solution u has components {density, x-momentum, y-momentum, (z-momentum), energy}.
   // These are stored contiguously in the BlockVector u_block.
   Array<int> offsets(num_equation + 1);
   for (int k = 0; k <= num_equation; k++) { offsets[k] = k * vfes.GetNDofs(); }
   BlockVector u_block(offsets);
   ParGridFunction sol(&vfes, u_block.GetData());

   ParGridFunction *saved_sol = NULL;

   // In case of restart read solution from the VisIt collection
   if (restart)
   {
      saved_sol = restart_dc.GetParField("restart_conservative_vars");
      if (!saved_sol)
      {
         if (myRank == 0)
         {
            cout << "The given data collection has no 'restart_conservative_vars' field."
                 << endl;
            exit(2);
         }
      }
      else
      {
         u_block.Update(*saved_sol,offsets);
         sol.MakeRef(&vfes, u_block, 0);
      }
   }
   else
   {
      // Initialize the state.
      VectorFunctionCoefficient u0(num_equation, InitialCondition);
      sol.ProjectCoefficient(u0);
      // Save the state to the VisIt
      restart_dc.RegisterField("restart_conservative_vars",&sol);
      restart_dc.Save();
   }

   if (myRank == 0) cout << "project sol OK\n";


   // 5. Set up the nonlinear form corresponding to the DG discretization of the
   //    flux divergence, and assemble the corresponding mass matrix.
   MixedBilinearForm Aflux(&dfes, &fes);
   Aflux.AddDomainIntegrator(new DomainIntegrator(dim));
   Aflux.Assemble();

   // cout << "after domain integr\n";

   ParNonlinearForm A(&vfes);
   RiemannSolver* rsolver = NULL;
   switch (riemann_solver_type)
   {
      case 0: rsolver = new RiemannSolverRusanov(); break;
      case 1: rsolver = new RiemannSolverLLF(); break;
      case 2: rsolver = new RiemannSolverHLL(); break;
      case 3: rsolver = new RiemannSolverHLLC(); break;
      default:
         cout << "Unknown Riemann solver type: " << riemann_solver_type << '\n';
         return 1;
   }
   A.AddInteriorFaceIntegrator(new FaceIntegrator(*rsolver, dim));

   Array<Array<int>> bdr_markers;
   SetBoundaryConditions(A, *pmesh, *rsolver, dim, bdr_markers);

   // cout << "after face interg\n";

   // 6. Define the time-dependent evolution operator describing the ODE
   //    right-hand side, and perform time-integration (looping over the time
   //    iterations, ti, with a time-step dt).
   FE_Evolution euler(vfes, A, Aflux.SpMat());

   // cout << "after fe evolution\n";




   // Determine the minimum element size.
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

   // 7. Define limiter and troubled cells indicator

   // 7.1. Define additional block vector for indicator values to have a possibility of visualisation
   DG_FECollection fec_const(0, dim);
   ParFiniteElementSpace fes_const(pmesh, &fec_const);
   Array<int> offsets_const(num_equation + 1);
   for (int k = 0; k <= num_equation; k++) { offsets_const[k] = k * fes_const.GetNDofs(); }
   BlockVector indicatorData(offsets_const);

   // 7.2. Define indicator object
   Indicator *ind = NULL;
   switch (indicator_type)
   {
      case 0: ind = new IndicatorNowhere(&vfes, offsets, dim, indicatorData); break;
      case 1: ind = new IndicatorEverywhere(&vfes, offsets, dim, indicatorData); break;
      case 2: ind = new IndicatorBJ(&vfes, offsets, dim, indicatorData); break;
      case 3: ind = new IndicatorShu(&vfes, offsets, dim, indicatorData); break;
      
      default:
         cout << "Unknown indicator type: " << indicator_type << '\n';
         return 1;
   }  

   // 7.3. Define limiter object
   Limiter *l = NULL;
   switch (limiter_type)
   {
      case 1: l = new LimiterFinDiff(*ind, &vfes, offsets, dim); break;
      case 2: l = new LimiterMultiplier(*ind, &vfes, offsets, dim); break;
      default:
         cout << "Unknown limiter type: " << limiter_type << '\n';
         return 1;
   }

   // 8. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.

   // int s = 2;
   // double* a = new double [3];
   // a[0] = 0.0;
   // a[1] = 1.0; 
   // a[2] = 0.0;
   // double* b = new double [2];
   // b[0] = 0.5;
   // b[1] = 0.5;
   // double* c = new double [2];
   // c[0] = 0.0;
   // c[1] = 1.0;

   int s = 3;
   double* a = new double [6];
   a[0] = 0.0;
   a[1] = 0.5; 
   a[2] = 0.0;
   a[3] = 0.0;
   a[4] = 1.0; 
   a[5] = 0.0;
   double* b = new double [3];
   b[0] = 0.1666666666666666;
   b[1] = 0.6666666666666666;
   b[2] = 0.1666666666666666;
   double* c = new double [2];
   c[0] = 0.0;
   c[1] = 0.5;
   c[2] = 1.0;
   ODESolver *ode_solver = new ExplicitRKLimitedSolver(s,a,b,c,*l);
   

   // 9. LIMIT INITIAL CONDITIONS
   // std::cout << "before initial limit " << sol[7052] << endl;
   pmesh->ExchangeFaceNbrData(); 
   sol.ExchangeFaceNbrData();
   l->update(sol);

   if (myRank == 0) cout << "limit sol OK\n";


   // 10. Create data collection for solution output: either VisItDataCollection for
   // ascii data files, or SidreDataCollection for binary data files.
   DataCollection *dc = NULL;

   ParGridFunction rhok(&fes, u_block.GetBlock(0));
   ParGridFunction mom(&dfes, u_block.GetData() + offsets[1]);
   ParGridFunction energy(&fes, u_block.GetBlock(dim+1));

   ParGridFunction rhoInd(&fes_const, indicatorData.GetBlock(0));
   ParGridFunction rhoUInd(&fes_const, indicatorData.GetBlock(1));
   ParGridFunction rhoVInd(&fes_const, indicatorData.GetBlock(2));
   ParGridFunction rhoWInd;
   if (dim == 3) rhoWInd = ParGridFunction(&fes_const, indicatorData.GetBlock(3));
   ParGridFunction EInd(&fes_const, indicatorData.GetBlock(dim+1));

   ParaViewDataCollection *pd = NULL;
   if (paraview)
   {
      pd = new ParaViewDataCollection("PV", pmesh);
      pd->RegisterField("mom", &mom);
      pd->RegisterField("rho", &rhok);
      pd->RegisterField("energy", &energy);
      pd->RegisterField("rhoInd", &rhoInd);
      pd->RegisterField("rhoUInd", &rhoUInd);
      pd->RegisterField("rhoVInd", &rhoVInd);
      if (dim == 3) pd->RegisterField("rhoWInd", &rhoWInd);
      pd->RegisterField("EInd", &EInd);

      pd->SetLevelsOfDetail(1);
      pd->SetCycle(0);
      pd->SetTime(0.0);
      pd->Save();
   }

   // Start the timer.
   tic_toc.Clear();
   tic_toc.Start();

   double t = restart ? restart_dc.GetTime() : 0.0;
   euler.SetTime(t);
   ode_solver->Init(euler);
   // cout << "after init euler\n";
   socketstream sout;

   if (!restart)
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


   // Integrate in time.
   if (myRank == 0) cout << "START TIME CYCLE" << endl;

   bool done = false;
   double t_real = 0.0;

   for (int ti = restart ? restart_dc.GetCycle() : 0 ; !done; )
   {
      t_real = tic_toc.RealTime();
      double dt_real = min(dt, t_final - t);

      ode_solver->Step(sol, t, dt_real);

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

      // for (int i = 0; i < num_equation; i++)
      // {
      //    sol_i.MakeRef(&fes, sol, offsets[i]);
      //    avgs_i.MakeRef(&fes_avg, avgs, offsets_avg[i]);
      //    sol_i.GetElementAverages(avgs_i);
      // }

      ti++;

      t_real = tic_toc.RealTime() - t_real;

      done = (t >= t_final - 1e-8*dt);

      if (done || ti % vis_steps == 0)
      {
         MPI_Barrier(pmesh->GetComm());

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
         
      }

      // done = (t >= t_final - 1e-8*dt);
      // if (done || ti % vis_steps == 0)
      // {
      //    cout << "time step: " << ti << "\tdt = " << dt << "\tPhys time: " << t << "\tReal time: " << t_real << endl;
      //    if (visualization)
      //    {
      //       sout << "solution\n" << mesh << mom << flush;
      //    }
      //    if (visit)
      //    {
      //       dc->SetCycle(ti);
      //       dc->SetTime(t);
      //       dc->Save();
      //    }


      // }
      if (myRank == 0)
      {
         cout << "Time step: " << ti << "\tdt: " << dt << "\tPhys time: " << t << " \tReal dt time: " << t_real << endl;
      }
   }

   tic_toc.Stop();
   if (myRank == 0) { cout << " done, " << tic_toc.RealTime() << "s." << endl; }


   // Free the used memory.
   delete pd;
   delete l;
   delete ode_solver;

   delete rsolver;

   delete[] c;
   delete[] b;
   delete[] a;


   MPI_Finalize();

   return 0;
}
