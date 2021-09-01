#include "rs_hll.hpp"
#include "physics.hpp"



double RiemannSolverHLL::Eval(const Vector &state1, const Vector &state2,
                           const Vector &nor, Vector &flux, bool debug)
{
   const int dim = nor.Size();

   StateIsPhysicalSay(state1, dim);
   StateIsPhysicalSay(state2, dim);

   Rotate(const_cast<Vector&>(state1), nor, dim);
   Rotate(const_cast<Vector&>(state2), nor, dim);

   ComputeFluxF(state1, dim, flux1);
   ComputeFluxF(state2, dim, flux2);

   ComputeToroCharSpeeds(state1, state2, lambdaF, dim);
   const double maxE = max(fabs(lambdaF[0]), fabs(lambdaF[dim+1]));

   if (lambdaF[0] >= 0)
      flux = flux1;

   else if (lambdaF[dim+1] <= 0)
      flux = flux2;

   else
      for (int i = 0; i < num_equation; ++i)
      {
         flux(i) = (lambdaF[dim+1] * flux1(i) - lambdaF[0] * flux2(i) + lambdaF[0] * lambdaF[dim+1] * (state2(i) - state1(i))) / (lambdaF[dim+1] - lambdaF[0]);
      }

   InverseRotate(flux, nor, dim);

   return maxE;
}
