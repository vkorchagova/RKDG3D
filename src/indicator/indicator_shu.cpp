#include "indicator_shu.hpp"

// bool pleaseWriteMe = false;

IndicatorShu::IndicatorShu
(
   Averager& _avgr, 
   ParFiniteElementSpace* _fes,
   ParFiniteElementSpace* _fes_const, 
   const Array<int>& _offsets, 
   int _d
) : Indicator(_avgr, _fes, _fes_const, _offsets, _d)  
{
   maxFabsPj.SetSize(num_equation);
   sumFabsDiffExtrap.SetSize(num_equation);
};


double IndicatorShu::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil,
   const DenseMatrix& elfun1_mat
)
{
   maxFabsPj = -DEFAULT_LARGE_NUMBER;
   sumFabsDiffExtrap = 0.0;

   Vector el_uMean(num_equation);
   Vector el_uMean_extrap(num_equation);

   // std::cout << "=== iCell = " << iCell << std::endl;

   averager.computeStencilExtrapAveragesVector(stencil);

   // if (iCell == 3052)
   // {
   //   std::cout << "avgs_extrap:\n";
   // find max fabs of average values among stencil
   for (int k : stencil->cell_num)
   {
      averager.readElementAverageByNumber(iCell, el_uMean);
      averager.readElementExtrapAverageByNumber(k, el_uMean_extrap);

      // find max fabs of average values among stencil
      for (int iSol = 0; iSol < num_equation; ++iSol)
      {
         // mI[iSol] = el_uMean[iSol] < mI[iSol] ? el_uMean[iSol] : mI[iSol];
         maxFabsPj[iSol] = fabs(el_uMean[iSol]) > maxFabsPj[iSol] ? fabs(el_uMean[iSol]) : maxFabsPj[iSol];

         sumFabsDiffExtrap[iSol] += fabs(el_uMean[iSol] - el_uMean_extrap[iSol]);

         // std::cout << "diff = " << el_uMean[iSol] << '-' <<el_uMean_extrap[iSol] << "; bool = " << (double)(k != iCell) << std::endl;
      }
   }
   // std::cout << "...\n";
   // }

   // if (iCell == 3052)
   // {
   //   std::cout << "maxFabsPj:";
   //   maxFabsPj.Print(std::cout);
   //   std::cout << "sumFabsDiffExtrap:";
   //   sumFabsDiffExtrap.Print(std::cout);
   // }

   // set indicator values to the external field

   double indVal = (fabs(maxFabsPj[0]) < DEFAULT_SMALL_EPSILON && fabs(sumFabsDiffExtrap[0]) < DEFAULT_SMALL_EPSILON) 
         ? 
         1.0 
         : 
         std::min(1.0, DEFAULT_SHU_CK * maxFabsPj[0] / (sumFabsDiffExtrap[0] + DEFAULT_SMALL_EPSILON));

   setValue(iCell, indVal);

   return indVal;


   // for (int iEq = 0; iEq < num_equation; ++iEq)
   //   values.GetBlock(iEq)[iCell] = (sumFabsDiffExtrap[iEq] / (maxFabsPj[iEq] + eps) > Ck) ? 0.0 : 1.0;

   // if (iCell == 1) exit(1);
};
