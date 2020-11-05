#include "rs_basic.hpp"
#include "physics.hpp"


// Implementation of class RiemannSolver
RiemannSolver::RiemannSolver() :
   flux1(num_equation),
   flux2(num_equation) { }


void RiemannSolver::Rotate(Vector& state, const Vector& nor, int dim)
{
   Vector mom(state.GetData() + 1, dim);
   double mx, my, mz;

   if (dim == 2)
   {
      mx =   nor[0] * mom[0] + nor[1] * mom[1];
      my = - nor[1] * mom[0] + nor[0] * mom[1];
      mom[0] = mx;
      mom[1] = my;
   }
   else if (dim == 3)
   {
      const double sinPsi = nor[2];
      const double cosPsi = sqrt(1 - sinPsi*sinPsi);
      const double sinTheta = fabs(cosPsi) > 1e-6 ? nor[1] / cosPsi : 1.0;
      const double cosTheta = fabs(cosPsi) > 1e-6 ? nor[0] / cosPsi : 0.0;

      mx =   cosPsi * cosTheta * mom[0] + cosPsi * sinTheta * mom[1] + sinPsi * mom[2];
      my =          - sinTheta * mom[0]          + cosTheta * mom[1];
      mz = - sinPsi * cosTheta * mom[0] - sinPsi * sinTheta * mom[1] + cosPsi * mom[2];

      mom[0] = mx;
      mom[1] = my;
      mom[2] = mz;
   } 
}

void RiemannSolver::InverseRotate(Vector& state, const Vector& nor, int dim)
{
   Vector mom(state.GetData() + 1, dim);
   double mx, my, mz;

   if (dim == 2)
   {
      mx = nor[0] * mom[0] - nor[1] * mom[1];
      my = nor[1] * mom[0] + nor[0] * mom[1];
      mom[0] = mx;
      mom[1] = my;
   }
   else if (dim == 3)
   {
      const double sinPsi = nor[2];
      const double cosPsi = sqrt(1 - sinPsi*sinPsi);
      const double sinTheta = fabs(cosPsi) > 1e-6 ? nor[1] / cosPsi : 1.0;
      const double cosTheta = fabs(cosPsi) > 1e-6 ? nor[0] / cosPsi : 0.0;

      mx = cosPsi * cosTheta * mom[0] - sinTheta * mom[1] - sinPsi * cosTheta * mom[2];
      my = cosPsi * sinTheta * mom[0] + cosTheta * mom[1] - sinPsi * sinTheta * mom[2];
      mz =            sinPsi * mom[0]                                + cosPsi * mom[2];

      mom[0] = mx;
      mom[1] = my;
      mom[2] = mz;
   } 
}