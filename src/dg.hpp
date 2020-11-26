#include "fe_evolution.hpp"
#include "domain_integrator.hpp"
#include "face_integrator.hpp"
#include "boundary_integrator_wall.hpp"
#include "boundary_integrator_open.hpp"
#include "boundary_integrator_constant.hpp"
#include "rs_rusanov.hpp"
#include "rs_hll.hpp"
#include "rs_hllc.hpp"
#include "rs_llf.hpp"
#include "rk_explicit_limited.hpp"
#include "limiter_findiff.hpp"
#include "limiter_multiplier.hpp"
#include "indicator_nowhere.hpp"
#include "indicator_everywhere.hpp"
#include "indicator_bj.hpp"
#include "indicator_shu.hpp"



// Initial condition
void InitialCondition(const Vector &x, Vector &y)
{
   const int dim = x.Size();

   double den = 0.0;
   double velX = 0.0;
   double velY = 0.0;
   double velZ = 0.0;
   double pres = 0.0;
   double energy = 0.0;
   double vel2 = 0.0;
   double shrinv1 = 1.0 / (specific_heat_ratio - 1.);

   if (problem == 1)
   { // Circle Sod problem
      
      double radius = 0.4;
      const double xc = 0.0, yc = 0.0, zc = 0.0;

      // double velX = 0.0;
      // double velY = 0.0;
      // double velZ = 0.0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;

         den  = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) < radius*radius ? 1.0 : 0.125;// * (x(0) + 1.0);
         pres = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) < radius*radius ? 1.0 : 0.1;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
      else if (dim == 3)
      {  
         vel2 = velX * velX + velY * velY + velZ * velZ;

         den  = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) + (x(2) - zc)*(x(2) - zc) < radius*radius ? 1.0 : 0.125;
         pres = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) + (x(2) - zc)*(x(2) - zc) < radius*radius ? 1.0 : 0.1;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
   }
   else if (problem == 2)
   { // Forward Step problem
      
      velX = 0.5;//0.675;//1.65;
      velY = 0.0;
      velZ = 0.0;

      den = 1.0;
      pres = 1.0 / specific_heat_ratio;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
      else if (dim == 3)
      {         
         mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
         vel2 = velX * velX + velY * velY + velZ * velZ;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
   }
   else if (problem == 3)
   { // Shu Osher problem
      
      const double xc = 0.125;

      
      velY = 0.0;
      velZ = 0.0;


      den  = x(0) < xc ? 3.857143 : 1 + 0.2*sin(8*2.0 * 3.14159265*x(0));// * (x(0) + 1.0);
      pres = x(0) < xc ? 10.33333 : 1;
      velX = x(0) < xc ? 2.629369 : 0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
      else if (dim == 3)
      {         
         mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
         vel2 = velX * velX + velY * velY + velZ * velZ;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
   }
   else if (problem == 4)
   { // Double Mach problem
      
      const double xc = 0.15;

      
      velY = 0.0;
      velZ = 0.0;


      den  = x(0) < xc ? 8.0 : 1.4;
      pres = x(0) < xc ? 116.518 : 1.0;
      velX = x(0) < xc ? 8.25 : 0.0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
      else if (dim == 3)
      {         
         mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
         vel2 = velX * velX + velY * velY + velZ * velZ;
         energy = shrinv1 * pres / den + 0.5 * vel2;
      }
   }
   else if (problem == 5)
   { // Sod Covolume problem

      covolume_constant = 0.001;
      specific_heat_ratio = 1.3;
      
      const double xc = 0.4;

      // double velX = 0.0;
      // double velY = 0.0;
      // double velZ = 0.0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;

         den  = x(0) < xc ? 100.0 : 1.0;// * (x(0) + 1.0);
         pres = x(0) < xc ? 1e8 : 1e5;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
      else if (dim == 3)
      {  
         vel2 = velX * velX + velY * velY + velZ * velZ;

         den  = x(0) < xc ? 100 : 1.0;
         pres = x(0) < xc ? 1e8 : 1e5;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
   }
   else if (problem == 6)
   { // State Contact problem
      
      const double xc = 0.5;

      double velX = 0.0;
      double velY = 0.0;
      double velZ = 0.0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;

         den  = x(0) < xc ? 1.4 : 1.0;// * (x(0) + 1.0);
         pres = 1.0;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
      else if (dim == 3)
      {  
         vel2 = velX * velX + velY * velY + velZ * velZ;

         den  = x(0) < xc ? 1.4 : 1.0;// * (x(0) + 1.0);
         pres = 1.0;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
   }
   else if (problem == 7)
   { // Strip Sod problem
      
      double stripLen = 1.0;
      const double xc = 0.0;

      // double velX = 0.0;
      // double velY = 0.0;
      // double velZ = 0.0;

      if (dim == 2)
      {
         vel2 = velX * velX + velY * velY;

         den  = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.125;// * (x(0) + 1.0);
         pres = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.1;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
      else if (dim == 3)
      {  
         vel2 = velX * velX + velY * velY + velZ * velZ;

         den  = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.125;
         pres = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.1;
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
   }
   else if (problem == 8)
   { // AstroTest

      specific_heat_ratio = 5.0/3.0;

      double u = 1.02; // radial velocity
      double kTilde = 5.5e6;
      double M0 = 1e30;
      double G0 = 6.67e-11;
      double Ror = 696e6;

      double k = kTilde * pow(M0, specific_heat_ratio - 2.0) / G0 / pow(Ror, 3.0*specific_heat_ratio - 4.0);
      double frackgamma = (specific_heat_ratio - 1.0) / k / specific_heat_ratio;

        // source = [=](const numvector<double, dimPh> sol, const Point& r)
        // {
        //     return numvector<double, 5> \
        //     { \
        //         0.0, \
        //         - x(0) / pow(x.Norml2(),3), \
        //         - x(1) / pow(x.Norml2(),3), \
        //         0.0, \
        //         0.0 \
        //     }; 
        // };

      if (dim == 2)
      {
         velX = - u * x(1) / pow( x.Norml2() * x.Norml2(), 0.75);
         velY = u * x(0) / pow( x.Norml2() * x.Norml2(), 0.75);

         vel2 = velX * velX + velY * velY;

         den  = pow( frackgamma * (u*u - 1.0) * (x.Norml2() - 1.0) / x.Norml2(), 1.0 / (specific_heat_ratio - 1.0));
         pres = k * pow( frackgamma * (u*u - 1.0) * (x.Norml2() - 1.0) / x.Norml2(), specific_heat_ratio / (specific_heat_ratio - 1.0));
         energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
      }
      else
      {
         cout << "wrong dim for AstroTest" << endl;
      }
   }
   else
   {
      mfem_error("Cannot recognize problem."
                 "Options are: 1 - Circle Sod problem, 2 - ForwardStep, 3 - Shu-Osher problem, 4 - DoubleMach problem, 5 - Covolume Sod problem, 6 - Contact State problem, 7 - StripSod");
   }   

   // set conservative variables
   if (dim == 2)
   {
      y(0) = den;
      y(1) = den * velX;
      y(2) = den * velY;
      y(3) = den * energy;
   }
   else if (dim == 3)
   {
      y(0) = den;
      y(1) = den * velX;
      y(2) = den * velY;
      y(3) = den * velZ;
      y(4) = den * energy;
   }
}

// Initial condition
void SetBoundaryConditions(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim, Array<Array<int>>& bdr_markers)
{
   // boundary attributes: active is 1, inactive is 0
   // Array<int> open_bdr_markers(mesh.bdr_attributes.Max());
   // Array<int> wall_bdr_markers(mesh.bdr_attributes.Max());
   // Array<int> inlet_bdr_markers(mesh.bdr_attributes.Max());

   // open_bdr_markers = 0; 
   // wall_bdr_markers = 0;
   // inlet_bdr_markers = 0;

   if (problem == 1 || problem == 5 || problem == 6 || problem == 7 || problem == 8) //CircleSod or Covolume or StateContact or StripSod
   {
      bdr_markers.SetSize(1);
      bdr_markers[0].SetSize(mesh.bdr_attributes.Max());
      bdr_markers[0][0] = 1; // let's set all boundaries opened
      bdr_markers[0][1] = 1;
      bdr_markers[0][2] = 1;
      bdr_markers[0][3] = 1;
      A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[0]);

   }
   else if (problem == 2) //ForwardStep
   {
      bdr_markers.SetSize(3); // three types of bcs
      for (int iBdr = 0; iBdr < 3; ++iBdr)
      {
         bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
         bdr_markers[iBdr] = 0;
      }

      bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
      bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
      bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (Obstacle)
      bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (Wall)

      double inletRho = 1.0;
      double inletU = 0.5;//0.675;//1.65; // opposite to normal of left boundary
      double inletV = 0.0;
      double inletP = 1.0 / specific_heat_ratio;

      Vector inletVars(dim + 2);

      if (dim == 3) mfem_error("Cannot initialize 3D problem for ForwardStep conditions");

      inletVars(0) = inletRho;
      inletVars(1) = inletRho * inletU;
      inletVars(2) = inletRho * inletV;
      inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


      A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[1]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
   }
   else if (problem == 3) //ShuOSher
   {
      bdr_markers.SetSize(3); // three types of bcs
      for (int iBdr = 0; iBdr < 3; ++iBdr)
      {
         bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
         bdr_markers[iBdr] = 0;
      }

      bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
      bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
      bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (Obstacle)
      bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (Wall)

      double inletRho = 3.857143;
      double inletU = 2.629369; // opposite to normal of left boundary
      double inletV = 0.0;
      double inletP = 10.33333;

      Vector inletVars(dim + 2);

      if (dim == 3) mfem_error("Cannot initialize 3D problem for ShuOsher conditions");

      inletVars(0) = inletRho;
      inletVars(1) = inletRho * inletU;
      inletVars(2) = inletRho * inletV;
      inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


      A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[1]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
   }
   else if (problem == 4) //DoubleMach
   {
      bdr_markers.SetSize(3); // three types of bcs
      for (int iBdr = 0; iBdr < 3; ++iBdr)
      {
         bdr_markers[iBdr].SetSize(mesh.bdr_attributes.Max());
         bdr_markers[iBdr] = 0;
      }

      bdr_markers[0][0] = 1; // boundary attribute 1 is Dirichlet (Inlet)
      bdr_markers[1][1] = 1; // boundary attribute 2 is Open      (Outlet)
      bdr_markers[2][2] = 1; // boundary attribute 3 is Wall      (top)
      bdr_markers[2][3] = 1; // boundary attribute 4 is Wall      (bottom)
      //bdr_markers[2][4] = 1; // boundary attribute 5 is Wall      (diag)

      double inletRho = 8.0;
      double inletU = 8.25; // opposite to normal of left boundary
      double inletV = 0.0;
      double inletP = 116.518;

      Vector inletVars(dim + 2);

      if (dim == 3) mfem_error("Cannot initialize 3D problem for ShuOsher conditions");

      inletVars(0) = inletRho;
      inletVars(1) = inletRho * inletU;
      inletVars(2) = inletRho * inletV;
      inletVars(3) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;


      A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[2]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[1]);
      A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[0]);
   }
   else
   {
      mfem_error("Cannot recognize problem to set boundary conditions."
                 "Options are: 1 - Circle Sod problem, 2 - ForwardStep, 3 - Shu-Osher problem, 4 - DoubleMach problem, 5 - Covolume Sod problem, 6 - StateContact, 7 - StripSod");
   }


}
