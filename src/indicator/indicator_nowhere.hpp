#ifndef INDICATOR_NOWHERE_H
#define INDICATOR_NOWHERE_H

#include "indicator.hpp"

using namespace std;
using namespace mfem;

extern const int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for checking discontinuities
/// for the DG slopes
///
class IndicatorNowhere : public Indicator
{   

public:

   /// Constructor
   IndicatorNowhere(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata)
      : Indicator(_fes, _offsets, _d, _idata) {};

   /// Destructor
   ~IndicatorNowhere() {};

   /// Find troubled cells
   virtual void checkDiscontinuity(const int iCell, 
      const Stencil* stencil, 
      const ParGridFunction* uMean, 
      const DenseMatrix& elfun1_mat,
      ParGridFunction &x) 
   { for (int iEq = 0; iEq < num_equation; ++iEq) values.GetBlock(iEq)[iCell] = 0.0; };
};

#endif // INDICATOR_NOWHERE_H

