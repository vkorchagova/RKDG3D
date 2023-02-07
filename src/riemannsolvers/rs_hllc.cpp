#include "rs_hllc.hpp"
#include "physics.hpp"

RiemannSolverHLLC::RiemannSolverHLLC() : RiemannSolver() 
{
   D.SetSize(num_equation);
   D = 0.0;
   FStar.SetSize(num_equation);
   FStar = 0.0;
   UStar.SetSize(num_equation);
   UStar = 0.0;
};

void RiemannSolverHLLC::getUStar (
   const Vector& state,
   const double pK,
   const double SK,
   const double cK,
   const double SStar,
   const int dim,
   Vector& UStar
) const
{
   double mult = state[0] * cK / (SK - SStar);

   double e = state[4] / state[0] + (SStar - state[1] / state[0]) * ( SStar + pK / state[0] / cK);

   UStar[0] = mult;
   UStar[1] = SStar * mult;

   for (int i = 2; i <=dim ; ++i)
      UStar[i] = state[i] * mult;

   UStar[dim+1] = e * mult;
}


void RiemannSolverHLLC::getFStar(
   const Vector& state, 
   const Vector& fK, 
   const double rhoL, 
   const double pK, 
   const double SK, 
   const double cK, 
   const double SStar,
   const int dim,
   Vector& FStar
)
{
   
   // var 0
   // getUStar(state, pK, SK, cK, SStar, dim, UStar);
   // for (int i = 0; i < state.Size(); ++i)
   //   FStar[i] = fK[i] + (UStar[i] - state[i]) * SK;

   

   // var 2
//   return ( SStar * (state*SK - fK) + D * pK * SK ) * (1.0 / (SK - SStar));

   // var 1
   D[1] = 1.0;
   D[dim+1] = SStar;
   for (int i = 0; i < state.Size(); ++i)
      FStar[i] = ( SStar * (state[i]*SK - fK[i]) + D[i] * SK * (pK + rhoL * cK * (SStar - state[1] / state[0])) ) * (1.0 / (SK - SStar));
}



double RiemannSolverHLLC::Eval(const Vector &state1, const Vector &state2,
                              const Vector &nor, Vector &flux, bool debug)
{
   const int dim = nor.Size();

   Rotate(const_cast<Vector&>(state1), nor, dim);
   Rotate(const_cast<Vector&>(state2), nor, dim);

   GetPrimitiveFromConservative(state1, primState1);
   GetPrimitiveFromConservative(state2, primState2);

   if (!StateIsPhysicalSay(state1, primState1, dim)) { std::cout << "Found in state 1 on proc #" << myRank << std::endl; return -1;};
   if (!StateIsPhysicalSay(state2, primState2, dim)) { std::cout << "Found in state 2 on proc #" << myRank << std::endl; return -1;};

   ComputeFluxF(state1, primState1, dim, flux1);
   ComputeFluxF(state2, primState2, dim, flux2);

   // if (debug)
   // {
   //   std::cout << "\tstate1 = ";
   //   state1.Print(std::cout);
   //   std::cout << "\tstate2 = ";
   //   state2.Print(std::cout);
   //   std::cout << "\tflux1 = ";
   //   flux1.Print(std::cout);
   //   std::cout << "\tflux2 = ";
   //   flux2.Print(std::cout);
   // }

   ComputeToroCharSpeeds(state1, state2, primState1, primState2, lambdaF, dim);
   const double maxE = std::max(fabs(lambdaF[0]), fabs(lambdaF[dim+1]));


   double SL = lambdaF[0];
   double SR = lambdaF[dim+1];

   // double SL = std::min(lambdaF[0],lambdaF[dim+1]);
   // double SR = std::max(lambdaF[0],lambdaF[dim+1]);

   // if (debug)
   // {
   //   std::cout << "SL = " << SL << std::endl;
   //   std::cout << "SR = " << SR << std::endl;
   // }

   if (SL >= 0)
      flux = flux1;

   else if (SR <= 0)
      flux = flux2;

   else
   {
      double pLeft = primState1[dim+1]; //ComputePressure(state1, dim);
      double pRight = primState2[dim+1]; // ComputePressure(state2, dim);

      double cLeft = SL - primState1[1]; //state1[1] / state1[0];
      double cRight = SR - primState2[1]; //state2[1] / state2[0];

      double sStar = (pRight - pLeft + state1[1] * cLeft - state2[1] * cRight) / \
                        (state1[0] * cLeft - state2[0] * cRight);
      // if (debug)
      // {
      //   std::cout << "pLeft = " << pLeft << std::endl;
      //   std::cout << "pRight = " << pRight << std::endl;
      //   std::cout << "cLeft = " << cLeft << std::endl;
      //   std::cout << "cRight = " << cRight << std::endl;
      //   std::cout << "sStar = " << sStar << std::endl;
      // }

      if (sStar > 0.0)
         getFStar( state1, flux1, state1[0],  pLeft, SL,  cLeft, sStar, dim, flux);

      else
         getFStar( state2, flux2, state2[0], pRight, SR, cRight, sStar, dim, flux);

      // if (debug)
      // {
      //   std::cout << "flux = ";
      //   flux.Print(std::cout);
      // }

   }

   InverseRotate(flux, nor, dim);
   // if (debug)
   // {
   //   std::cout << "flux after rot = ";
   //   flux.Print(std::cout);
   // }


   return maxE;
}
