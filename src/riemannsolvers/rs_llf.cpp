#include "rs_llf.hpp"
#include "physics.hpp"



double RiemannSolverLLF::Eval(const Vector &state1, const Vector &state2,
                           const Vector &nor, Vector &flux, bool debug)
{
   // NOTE: nor in general is not a unit normal
   const int dim = nor.Size();

   // state1.Print(cout);
   // state2.Print(cout);
   StateIsPhysicalSay(state1, dim);
   StateIsPhysicalSay(state2, dim);

   
   Rotate(const_cast<Vector&>(state1), nor, dim);
   Rotate(const_cast<Vector&>(state2), nor, dim);

   ComputeFluxF(state1, dim, flux1);
   ComputeFluxF(state2, dim, flux2);

   ComputeEinfeldtCharSpeeds(state1, state2, lambdaF, dim);

   const double maxE = max(fabs(lambdaF[0]), fabs(lambdaF[dim+1]));

   for (int i = 0; i < num_equation; ++i)
   {
      flux(i) = 0.5 * (flux1(i) + flux2(i))
                - 0.5 * maxE * (state2(i) - state1(i));
   }

   InverseRotate(flux, nor, dim);

   return maxE;
}
