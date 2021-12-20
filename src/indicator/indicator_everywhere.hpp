#ifndef INDICATOR_EVERYWHERE_H
#define INDICATOR_EVERYWHERE_H

#include "indicator.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for checking discontinuities
/// for the DG slopes
///
class IndicatorEverywhere : public Indicator
{   

public:

   /// Constructor
   IndicatorEverywhere(Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata)
      : Indicator(_avgr, _fes,_offsets,_d, _idata) {};

   /// Destructor
   ~IndicatorEverywhere() {};

   /// Find troubled cells
   virtual void checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) 
   { for (int iEq = 0; iEq < num_equation; ++iEq) {values.GetBlock(iEq)[iCell] = 0.0; minValues[iCell] = 0.0;} };
};

#endif // INDICATOR_EVERYWHERE_H

