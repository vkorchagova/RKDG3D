#include "limiter_findiff.hpp"

void LimiterFinDiff::limit(const int iCell, const double ind_value, const double nDofs, DenseMatrix& elfun1_mat) 
{
   averager.readElementAverageByNumber(iCell, el_uMean);

   // replace solution to mean values
   for (int iEq = 0; iEq < num_equation; ++iEq)
      for (int iDof = 0; iDof < nDofs; ++iDof)
      {
        elfun1_mat(iDof, iEq) = ind_value < DEFAULT_INDICATOR_CORRECTION_THRESHOLD_VALUE ? el_uMean(iEq) : elfun1_mat(iDof, iEq);
      }
};   
