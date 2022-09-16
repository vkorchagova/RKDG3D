#ifndef INDICATOR_NOWHERE_H
#define INDICATOR_NOWHERE_H

#include "indicator.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Set all cells should NOT be be limited
///
class IndicatorNowhere : public Indicator
{   

public:

   /// Constructor
   IndicatorNowhere
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   ) : Indicator(_avgr, _fes, _fes_const, _offsets, _d) {};

   /// Destructor
   ~IndicatorNowhere() {};

   /// Find troubled cells
   virtual double checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil,
      const DenseMatrix& elfun1_mat
   ) 
   { 
      setValue(iCell, 1.0); return 1.0;
   };
};

#endif // INDICATOR_NOWHERE_H

