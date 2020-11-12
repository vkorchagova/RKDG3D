#include "indicator_shu.hpp"

// bool pleaseWriteMe = false;

IndicatorShu::IndicatorShu(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata) 
    : Indicator(_fes, _offsets, _d, _idata) 
{
   maxFabsPj.SetSize(num_equation);
   sumFabsDiffExtrap.SetSize(num_equation);

   // prepare place for the average values
   fec_avg_extrap = new DG_FECollection(0, dim);
   fes_avg_extrap = new ParFiniteElementSpace(mesh, fec_avg_extrap, num_equation);
   offsets_avg_extrap.SetSize(num_equation + 1);

   fes_avg_extrap_component = new ParFiniteElementSpace(mesh, fec_avg_extrap);

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg_extrap[k] = k * fes_avg_extrap->GetNDofs();

   u_block_avg_extrap = new BlockVector(offsets_avg_extrap);
   avgs_extrap = new ParGridFunction(fes_avg_extrap, u_block_avg_extrap->GetData());
};


void IndicatorShu::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil, 
   const ParGridFunction* uMean, 
   const DenseMatrix& elfun1_mat,
   ParGridFunction &x)
{
   maxFabsPj = -1e9;
   sumFabsDiffExtrap = 0.0;

   Vector el_uMean(num_equation);
   Vector el_uMean_extrap(num_equation);

   // cout << "=== iCell = " << iCell << endl;

   computeStencilExtrapAveragesVector(x, stencil, fes, fes_avg_extrap_component, fec_avg_extrap, offsets, offsets_avg_extrap, mesh, avgs_extrap);

   // if (iCell == 3052)
   // {
   //    cout << "avgs_extrap:\n";
   // find max fabs of average values among stencil
   for (int k : stencil->cell_num)
   {
      readElementAverageByNumber(iCell, mesh, uMean, el_uMean);
      readElementAverageByNumber(k, mesh, avgs_extrap, el_uMean_extrap);

      // find max fabs of average values among stencil
      for (int iSol = 0; iSol < num_equation; ++iSol)
      {
         // mI[iSol] = el_uMean[iSol] < mI[iSol] ? el_uMean[iSol] : mI[iSol];
         maxFabsPj[iSol] = fabs(el_uMean[iSol]) > maxFabsPj[iSol] ? fabs(el_uMean[iSol]) : maxFabsPj[iSol];

         sumFabsDiffExtrap[iSol] += fabs(el_uMean[iSol] - el_uMean_extrap[iSol]) * (double)(k != iCell);

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
      values.GetBlock(iEq)[iCell] = min(1.0, Ck * maxFabsPj[iEq] / (sumFabsDiffExtrap[iEq] + eps));

   // for (int iEq = 0; iEq < num_equation; ++iEq)
   //    values.GetBlock(iEq)[iCell] = (sumFabsDiffExtrap[iEq] / (maxFabsPj[iEq] + eps) > Ck) ? 0.0 : 1.0;

   // if (iCell == 1) exit(1);
};
