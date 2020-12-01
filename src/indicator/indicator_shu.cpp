#include "indicator_shu.hpp"

// bool pleaseWriteMe = false;

IndicatorShu::IndicatorShu(Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata) 
    : Indicator(_avgr, _fes, _offsets, _d, _idata) 
{
   maxFabsPj.SetSize(num_equation);
   sumFabsDiffExtrap.SetSize(num_equation);
};


void IndicatorShu::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil,
   const DenseMatrix& elfun1_mat
)
{
   maxFabsPj = -1e9;
   sumFabsDiffExtrap = 0.0;

   Vector el_uMean(num_equation);
   Vector el_uMean_extrap(num_equation);

   // cout << "=== iCell = " << iCell << endl;

   averager.computeStencilExtrapAveragesVector(stencil);

   // if (iCell == 3052)
   // {
   //    cout << "avgs_extrap:\n";
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

         // cout << "diff = " << el_uMean[iSol] << '-' <<el_uMean_extrap[iSol] << "; bool = " << (double)(k != iCell) << endl;
      }
   }
   // cout << "...\n";
   // }

   // if (iCell == 3052)
   // {
   //    cout << "maxFabsPj:";
   //    maxFabsPj.Print(cout);
   //    cout << "sumFabsDiffExtrap:";
   //    sumFabsDiffExtrap.Print(cout);
   // }

   // set indicator values to the external field

   for (int iEq = 0; iEq < num_equation; ++iEq)
      values.GetBlock(iEq)[iCell] = 
         (fabs(maxFabsPj[iEq]) < eps && fabs(sumFabsDiffExtrap[iEq]) < eps) 
         ? 
         1.0 
         : 
         min(1.0, Ck * maxFabsPj[iEq] / (sumFabsDiffExtrap[iEq] + eps));

   // for (int iEq = 0; iEq < num_equation; ++iEq)
   //    values.GetBlock(iEq)[iCell] = (sumFabsDiffExtrap[iEq] / (maxFabsPj[iEq] + eps) > Ck) ? 0.0 : 1.0;

   // if (iCell == 1) exit(1);
};
