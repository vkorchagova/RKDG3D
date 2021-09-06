#include "rs_rusanov.hpp"
#include "physics.hpp"



double RiemannSolverRusanov::Eval(const Vector &state1, const Vector &state2,
                           const Vector &nor, Vector &flux, bool debug)
{
   // NOTE: nor in general is not a unit normal
   const int dim = nor.Size();

   if (!StateIsPhysicalSay(state1, dim)) exit(1);
   if (!StateIsPhysicalSay(state2, dim)) exit(1);

   const double maxE1 = ComputeMaxCharSpeed(state1, dim);
   const double maxE2 = ComputeMaxCharSpeed(state2, dim);

   // state1.Print(cout);
   // state2.Print(cout);

   const double maxE = max(maxE1, maxE2);

   ComputeFluxDotN(state1, nor, flux1);
   ComputeFluxDotN(state2, nor, flux2);

   double normag = 0;
   for (int i = 0; i < dim; i++)
   {
      normag += nor(i) * nor(i);
   }
   normag = sqrt(normag);

   for (int i = 0; i < num_equation; i++)
   {
      flux(i) = 0.5 * (flux1(i) + flux2(i)) / normag
                - 0.5 * maxE * (state2(i) - state1(i));
   }

   return maxE;
}
