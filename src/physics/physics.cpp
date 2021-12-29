
#include "physics.hpp"


// Check that the state is physical
bool StateIsPhysical(const Vector &state, const int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   if (den < 0 || den != den)
   {
      // cout << "Negative density: ";
      // for (int i = 0; i < state.Size(); i++)
      // {
      //    cout << state(i) << " ";
      // }
      // cout << endl;
      return false;
   }
   if (den_energy <= 0 || den_energy != den_energy)
   {
      // cout << "Negative energy: ";
      // for (int i = 0; i < state.Size(); i++)
      // {
      //    cout << state(i) << " ";
      // }
      // cout << endl;
      return false;
   }

   double den_vel2 = 0;
   for (int i = 0; i < dim; i++) { den_vel2 += den_vel(i) * den_vel(i); }
   den_vel2 /= den;

   const double pres = ComputePressure(state, dim);

   if (pres <= 0 || pres != pres)
   {
      // cout << "Negative pressure: " << pres << ", state: ";
      // for (int i = 0; i < state.Size(); i++)
      // {
      //    cout << state(i) << " ";
      // }
      // cout << endl;
      return false;
   }
   return true;
}

// Check that the state is physical
bool StateIsPhysicalSay(const Vector &state, const int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   if (den < 0 || den != den)
   {
      cout << "Rank #" << myRank << ": " << "Negative density: ";
      for (int i = 0; i < state.Size(); i++)
      {
         cout << state(i) << " ";
      }
      cout << endl;
      return false;
   }
   if (den_energy <= 0 || den_energy != den_energy)
   {
      cout << "Rank #" << myRank << ": " << "Negative energy: ";
      for (int i = 0; i < state.Size(); i++)
      {
         cout << state(i) << " ";
      }
      cout << endl;
      return false;
   }

   double den_vel2 = 0;
   for (int i = 0; i < dim; i++) { den_vel2 += den_vel(i) * den_vel(i); }
   den_vel2 /= den;

   const double pres = (specific_heat_ratio - 1.0) * (den_energy - 0.5 * den_vel2);

   if (pres <= 0 || pres != pres)
   {
      cout << "Rank #" << myRank << ": " << "Negative pressure: " << pres << ", state: ";
      for (int i = 0; i < state.Size(); i++)
      {
         cout << state(i) << " ";
      }
      cout << endl;
      return false;
   }
   return true;
}


// Pressure (EOS) computation
double ComputePressure(const Vector &state, int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   double den_vel2 = 0;
   for (int d = 0; d < dim; d++) { den_vel2 += den_vel(d) * den_vel(d); }
   den_vel2 /= den;

   return (specific_heat_ratio - 1.0) * (den_energy - 0.5 * den_vel2) / (1.0 - den * covolume_constant);
}

// Pressure (EOS) computation
double ComputeTemperature(const Vector &state, int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   double den_vel2 = 0;
   for (int d = 0; d < dim; d++) { den_vel2 += den_vel(d) * den_vel(d); }
   den_vel2 /= den;

   return (specific_heat_ratio - 1.0) * (den_energy - 0.5 * den_vel2) / den / gas_constant;
}

double ComputeEnergy(double rho, double u, double v, double w, double p)
{
   return p / (specific_heat_ratio - 1.0) * (1.0 - rho * covolume_constant) + 0.5 * rho * (u*u + v*v + w*w);
}

// Sound speed (EOS) computation
double ComputeSoundSpeed(const Vector &state, int dim)
{
   const double den = state(0);
   const double p = ComputePressure(state, dim);

   return sqrt( specific_heat_ratio * p / den / (1.0 - den * covolume_constant));
}

// Compute the vector flux F(u)
void ComputeFlux(const Vector &state, int dim, DenseMatrix &flux)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   MFEM_ASSERT(StateIsPhysical(state, dim), "");

   const double pres = ComputePressure(state, dim);

   for (int d = 0; d < dim; d++)
   {
      flux(0, d) = den_vel(d);
      for (int i = 0; i < dim; i++)
      {
         flux(1+i, d) = den_vel(i) * den_vel(d) / den;
      }
      flux(1+d, d) += pres;
   }

   const double H = (den_energy + pres) / den;
   for (int d = 0; d < dim; d++)
   {
      flux(1+dim, d) = den_vel(d) * H;
   }
}

// Compute the scalar F(u).n
void ComputeFluxDotN(const Vector &state, const Vector &nor,
                     Vector &fluxN)
{
   // NOTE: nor in general is not a unit normal
   const int dim = nor.Size();
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);
   const double den_energy = state(1 + dim);

   MFEM_ASSERT(StateIsPhysical(state, dim), "");

   const double pres = ComputePressure(state, dim);

   double den_velN = 0;
   for (int d = 0; d < dim; d++) { den_velN += den_vel(d) * nor(d); }

   fluxN(0) = den_velN;
   for (int d = 0; d < dim; d++)
   {
      fluxN(1+d) = den_velN * den_vel(d) / den + pres * nor(d);
   }

   const double H = (den_energy + pres) / den;
   fluxN(1 + dim) = den_velN * H;
}

// Compute the scalar F(u).n
void ComputeFluxF(const Vector &state, const int dim,
                     Vector &flux)
{
   const double pres = ComputePressure(state, dim);

   flux(0) = state(1); //rho u
   flux(1) = state(1) * state(1) / state(0) + pres; // rho u^2 + p
   for (int d = 1; d < dim; d++)
      flux(d+1) = state(1) * state(d+1) / state(0); // rho u v, rho u w

   const double H = (state(1 + dim) + pres) / state(0);
   flux(1 + dim) = state(1) * H;



   // const double den = state(0);
   // const Vector den_vel(state.GetData() + 1, dim);
   // const double den_energy = state(1 + dim);

   // MFEM_ASSERT(StateIsPhysical(state, dim), "");

   

   // double den_velN = 0;
   // for (int d = 0; d < dim; d++) { den_velN += den_vel(d) * nor(d); }

   // fluxN(0) = den_velN;
   // for (int d = 0; d < dim; d++)
   // {
   //    fluxN(1+d) = den_velN * den_vel(d) / den + pres * nor(d);
   // }

   
}


// Compute the Mach number.
double ComputeM(const Vector &state, const int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);

   double den_vel2 = 0;
   for (int d = 0; d < dim; d++) { den_vel2 += den_vel(d) * den_vel(d); }
   den_vel2 /= den;

   // const double pres = ComputePressure(state, dim);
   const double sound = ComputeSoundSpeed(state, dim);
   const double vel = sqrt(den_vel2 / den);

   return vel/sound;
}

// Compute the maximum characteristic speed.
double ComputeMaxCharSpeed(const Vector &state, const int dim)
{
   const double den = state(0);
   const Vector den_vel(state.GetData() + 1, dim);

   double den_vel2 = 0;
   for (int d = 0; d < dim; d++) { den_vel2 += den_vel(d) * den_vel(d); }
   den_vel2 /= den;

   // const double pres = ComputePressure(state, dim);
   const double sound = ComputeSoundSpeed(state, dim);
   const double vel = sqrt(den_vel2 / den);

   return vel + sound;
}

// Compute Einfeldt averaged char speeds via two states
void ComputeEinfeldtCharSpeeds(const Vector &state1, const Vector &state2, Vector& lambdaF, const int dim)
{
   double cLeft = ComputeSoundSpeed(state1, dim);
   double cRight = ComputeSoundSpeed(state2, dim);

   double rhouLeft = state1[1];
   double rhouRight = state2[1];

   double uLeft = rhouLeft / state1[0];
   double uRight = rhouRight / state2[0];

   double sqrtRhoLeft = sqrt(state1[0]);
   double sqrtRhoRight = sqrt(state2[0]);
   double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

   double eta2 = 0.5 * sqrtRhoLeft * sqrtRhoRight / sqr(sqrtRhoLeft + sqrtRhoRight);
   double u_av = (rhouLeft / sqrtRhoLeft + rhouRight / sqrtRhoRight) / sumSqrtRho;
   double c_av = sqrt((sqrtRhoLeft * sqr(cLeft) + sqrtRhoRight * sqr(cRight)) / sumSqrtRho + \
      eta2 * sqr(uRight - uLeft));

   // Vector res(dim+2);
   lambdaF[0] = u_av - c_av;
   lambdaF[dim+1] = u_av + c_av;

   for (int i = 1; i < dim + 1; ++i)
      lambdaF[i] = u_av;

}

void TransformConservativeToPrimitive(Vector& state)
{
   const int dim = state.Size() - 2;

   for (int i = 1; i < dim+1; ++i)
      state[i] /= state[0];

   state[dim+1] = ComputePressure(state, dim);
}

// Compute Toro averaged char speeds via two states
void ComputeToroCharSpeeds(const Vector &state1, const Vector &state2, Vector& lambdaF, const int dim)
{
   double cLeft = ComputeSoundSpeed(state1, dim);
   double cRight = ComputeSoundSpeed(state2, dim);

   double rhouLeft = state1[1];
   double rhouRight = state2[1];

   double uLeft = rhouLeft / state1[0];
   double uRight = rhouRight / state2[0];

   double pLeft = ComputePressure(state1, dim);
   double pRight = ComputePressure(state2, dim);

   //  double pvrs = 0.5 * (pLeft + pRight) + \
   //    0.125 * (uLeft - uRight) * (state1[0] + state2[0]) * (cLeft + cRight);

   // double pStar = max(0.0, pvrs);

   double z = 0.5 * (specific_heat_ratio - 1.0) / specific_heat_ratio;

   double pStar = pow( 
      (cLeft + cRight - 0.5 * (specific_heat_ratio - 1.0) * (uRight - uLeft)) / \
      (cLeft / pow(pLeft, z) + cRight / pow(pRight, z)),
      1.0 / z);

   double qLeft = pStar > pLeft ? \
      sqrt(1.0 + 0.5 * (specific_heat_ratio + 1.0) * (pStar / pLeft - 1.0) / specific_heat_ratio) : \
      1.0;

   double qRight = pStar > pRight ? \
      sqrt(1.0 + 0.5 * (specific_heat_ratio + 1.0) * (pStar / pRight - 1.0) / specific_heat_ratio) : \
      1.0;

   // Vector res(dim+2);
   lambdaF[0] = uLeft - qLeft * cLeft;
   lambdaF[dim+1] = uRight + qRight * cRight;

   for (int i = 1; i < dim + 1; ++i)
      lambdaF[i] = uLeft;

}


// void ComputeU(const GridFunction& sol, const FiniteElementSpace& vfes, GridFunction& U)
// {
//     for (int i = 0; i < vfes.GetNE(); i++)
//     {
//         Vector zval;
//         Array<int> vdofs;
//         vfes.GetElementVDofs(i, vdofs);
//         z.GetSubVector(vdofs, zval);
//         // if (i == 50 || i == 1562)
//         // {
//         //     zval.Print(std::cout << "zval[" << i << "] after fluxes = " << endl);
//         // }
//     }
// }
