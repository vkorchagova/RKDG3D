#include "fe_evolution.hpp"
#include "domain_integrator.hpp"
#include "face_integrator.hpp"




// // Initial condition
// void SetBoundaryConditions(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim, Array<Array<int>>& bdr_markers)
// {
//    // boundary attributes: active is 1, inactive is 0
//    // Array<int> open_bdr_markers(mesh.bdr_attributes.Max());
//    // Array<int> wall_bdr_markers(mesh.bdr_attributes.Max());
//    // Array<int> inlet_bdr_markers(mesh.bdr_attributes.Max());

//    // open_bdr_markers = 0; 
//    // wall_bdr_markers = 0;
//    // inlet_bdr_markers = 0;

//    if (problem == 1 || problem == 5 || problem == 6 || problem == 7 || problem == 8) //CircleSod or Covolume or StateContact or StripSod
//    {
//       bdr_markers.SetSize(1);
//       bdr_markers[0].SetSize(mesh.bdr_attributes.Max());
//       bdr_markers[0][0] = 1; // let's set all boundaries opened
//       bdr_markers[0][1] = 1;
//       bdr_markers[0][2] = 1;
//       bdr_markers[0][3] = 1;
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[0]);

//    }
//    else if (problem == 2) //ForwardStep
//    {
//       bdr_markers.SetSize(3); // three types of bcs
//       for (int iBdr = 0; iBdr < 3; ++iBdr)
//       {
//          bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
//          bdr_markers[iBdr] = 0;
//       }

//       bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
//       bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
//       bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (Obstacle)
//       bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (Wall)

//       double inletRho = 1.0;
//       double inletU = -0.1;//0.675;//1.65; // opposite to normal of left boundary
//       double inletV = 0.0;
//       double inletP = 1.0 / specific_heat_ratio;

//       Vector inletVars(dim + 2);

//       if (dim == 3) mfem_error("Cannot initialize 3D problem for ForwardStep conditions");

//       inletVars(0) = inletRho;
//       inletVars(1) = inletRho * inletU;
//       inletVars(2) = inletRho * inletV;
//       inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[1]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
//    }
//    else if (problem == 3) //ShuOSher
//    {
//       bdr_markers.SetSize(3); // three types of bcs
//       for (int iBdr = 0; iBdr < 3; ++iBdr)
//       {
//          bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
//          bdr_markers[iBdr] = 0;
//       }

//       bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
//       bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
//       bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (Obstacle)
//       bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (Wall)

//       double inletRho = 3.857143;
//       double inletU = 2.629369; // opposite to normal of left boundary
//       double inletV = 0.0;
//       double inletP = 10.33333;

//       Vector inletVars(dim + 2);

//       if (dim == 3) mfem_error("Cannot initialize 3D problem for ShuOsher conditions");

//       inletVars(0) = inletRho;
//       inletVars(1) = inletRho * inletU;
//       inletVars(2) = inletRho * inletV;
//       inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[1]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
//    }
//    else if (problem == 4) //DoubleMach
//    {
//       bdr_markers.SetSize(3); // three types of bcs
//       for (int iBdr = 0; iBdr < 3; ++iBdr)
//       {
//          bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
//          bdr_markers[iBdr] = 0;
//       }

//       bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
//       bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
//       bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (top)
//       bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (bottom)
//       //bdr_markers[2][4] = 1; // boundary attribute 5 is Wall      (diag)

//       double inletRho = 8.0;
//       double inletU = 8.25; // opposite to normal of left boundary
//       double inletV = 0.0;
//       double inletP = 116.518;

//       Vector inletVars(dim + 2);

//       if (dim == 3) mfem_error("Cannot initialize 3D problem for ShuOsher conditions");

//       inletVars(0) = inletRho;
//       inletVars(1) = inletRho * inletU;
//       inletVars(2) = inletRho * inletV;
//       inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[1]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
//    }
//    else if (problem == 9) //Wing
//    {
//       bdr_markers.SetSize(3); // three types of bcs
//       for (int iBdr = 0; iBdr < 3; ++iBdr)
//       {
//          bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
//          bdr_markers[iBdr] = 0;
//       }

//       bdr_markers[1][0] = 1; // boundary attribute 1 is Open   (topAndBottom)
//       bdr_markers[0][1] = 1; // boundary attribute 2 is Inlet  (inlet)
//       bdr_markers[1][2] = 1; // boundary attribute 3 is Open   (outlet)
//       bdr_markers[2][3] = 1; // boundary attribute 4 is Wall   (Wing)

//       double inletRho = 1.176413;
//       double inletU = 577.2273;//0.675;//1.65; // opposite to normal of left boundary
//       double inletV = 0.0;
//       double inletP = 101325;

//       Vector inletVars(dim + 2);

//       if (dim == 3) mfem_error("Cannot initialize 3D problem for ForwardStep conditions");

//       inletVars(0) = inletRho;
//       inletVars(1) = inletRho * inletU;
//       inletVars(2) = inletRho * inletV;
//       inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


//       A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[1]);
//       A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
//    }
//    else
//    {
//       mfem_error("Cannot recognize problem to set boundary conditions."
//                  "Options are: 1 - Circle Sod problem, 2 - ForwardStep, 3 - Shu-Osher problem, 4 - DoubleMach problem, 5 - Covolume Sod problem, 6 - StateContact, 7 - StripSod, 8 - AstroTest, 9 - Wing");
//    }


// }
