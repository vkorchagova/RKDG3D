#include "rs_rusanov.hpp"
#include "physics.hpp"



double RiemannSolverRusanov::Eval(const Vector &state1, const Vector &state2,
                              const Vector &nor, Vector &flux, bool debug)
{
   // NOTE: nor in general is not a unit normal
   const int dim = nor.Size();

   GetPrimitiveFromConservative(state1, primState1);
   GetPrimitiveFromConservative(state2, primState2);

   if (!StateIsPhysicalSay(state1, primState1, dim)) { std::cout << "Found in state 1 on proc #" << myRank; return -1;};
   if (!StateIsPhysicalSay(state2, primState2, dim)) { std::cout << "Found in state 2 on proc #" << myRank; return -1;};
   
   const double maxE1 = ComputeMaxCharSpeed(state1, dim);
   const double maxE2 = ComputeMaxCharSpeed(state2, dim);

   // state1.Print(std::cout);
   // state2.Print(std::cout);

   const double maxE = std::max(maxE1, maxE2);

   ComputeFluxDotN(state1, nor, flux1);
   ComputeFluxDotN(state2, nor, flux2);

   double normag = 0;
   for (int i = 0; i < dim; ++i)
   {
      normag += nor(i) * nor(i);
   }
   normag = sqrt(normag);

   for (int i = 0; i < num_equation; ++i)
   {
      flux(i) = 0.5 * (flux1(i) + flux2(i)) / normag
                - 0.5 * maxE * (state2(i) - state1(i));
   }

   return maxE;
}
