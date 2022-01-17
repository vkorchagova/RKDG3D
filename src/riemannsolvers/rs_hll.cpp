#include "rs_hll.hpp"
#include "physics.hpp"



double RiemannSolverHLL::Eval(const Vector &state1, const Vector &state2,
                           const Vector &nor, Vector &flux, bool debug)
{
   const int dim = nor.Size();

   if (!StateIsPhysicalSay(state1, dim)) { cout << "Found in state 1 on proc #" << myRank; return -1;};
   if (!StateIsPhysicalSay(state2, dim)) { cout << "Found in state 2 on proc #" << myRank; return -1;};

   Rotate(const_cast<Vector&>(state1), nor, dim);
   Rotate(const_cast<Vector&>(state2), nor, dim);

   GetPrimitiveFromConservative(state1, primState1);
   GetPrimitiveFromConservative(state2, primState2);

   ComputeFluxF(state1, primState1, dim, flux1);
   ComputeFluxF(state2, primState2, dim, flux2);
   
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

   // if (debug)
   // {
      
   //    state1.Print(cout << std::setprecision(30) << " * state 1 = ");
   //    state2.Print(cout << std::setprecision(30) << " * state 2 = ");
   //    flux1.Print(cout << std::setprecision(30) << " * flux 1 = ");
   //    flux2.Print(cout << std::setprecision(30) << " * flux 2 = ");
   //    lambdaF.Print(cout << std::setprecision(30) << " * lambdaF = ");
   //    // cout << lambdaF[dim+1] << ' ' << flux1(3) << ' ' << lambdaF[0] << ' ' << flux2(3) << ' ' << lambdaF[0] << ' ' << lambdaF[dim+1] << ' ' << (state2(3) - state1(3)) << ' ' << (lambdaF[dim+1] - lambdaF[0]) << endl;
   //    cout << std::setprecision(30) << state2(3) << ' ' << state1(3) << ' ' << (state2(3) - state1(3))<< endl;
   // }

   return maxE;
}
