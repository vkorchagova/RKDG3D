#include "rs_llf.hpp"
#include "physics.hpp"

double RiemannSolverLLF::Eval(const Vector &state1, const Vector &state2,
                              const Vector &nor, Vector &flux, bool debug)
{
   // NOTE: nor in general is not a unit normal
   const int dim = nor.Size();

   if (debug)
   {
      state1.Print(std::cout << std::setprecision(18) << "state1 = ");
      state2.Print(std::cout << std::setprecision(18) << "state2 = ");
      nor.Print(std::cout << std::setprecision(18) << "nor = ");
   }

   Rotate(const_cast<Vector&>(state1), nor, dim);
   Rotate(const_cast<Vector&>(state2), nor, dim);

   if (debug)
   {
      state1.Print(std::cout << std::setprecision(18) << "state1 after rot = ");
      state2.Print(std::cout << std::setprecision(18) << "state2 after rot = ");
   }

   GetPrimitiveFromConservative(state1, primState1);
   GetPrimitiveFromConservative(state2, primState2);

   if (!StateIsPhysicalSay(state1, primState1, dim)) { std::cout << "Found in state 1 on proc #" << myRank << std::endl; return -1;};
   if (!StateIsPhysicalSay(state2, primState2, dim)) { std::cout << "Found in state 2 on proc #" << myRank << std::endl; return -1;};

   ComputeFluxF(state1, primState1, dim, flux1);
   ComputeFluxF(state2, primState2, dim, flux2);

   ComputeToroCharSpeeds(state1, state2, primState1, primState2, lambdaF, dim);

   const double maxE = std::max(fabs(lambdaF[0]), fabs(lambdaF[dim+1]));

   if (debug)
   {
      flux1.Print(std::cout << std::setprecision(18) << "flux1 = ");
      flux2.Print(std::cout << std::setprecision(18) << "flux2 = ");
      lambdaF.Print(std::cout << std::setprecision(18) << "lambdaF = ");
      std::cout << "maxE = " << maxE << std::endl;
   }


   for (int i = 0; i < num_equation; ++i)
   {
       flux(i) = 0.5 * (flux1(i) + flux2(i))
                 - 0.5 * maxE * (state2(i) - state1(i));
   }


   InverseRotate(flux, nor, dim);

   return maxE;
}

