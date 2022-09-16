#ifndef INDICATOR_EVERYWHERE_H
#define INDICATOR_EVERYWHERE_H

#include "indicator.hpp"


using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Set all cells should be limited
///
class IndicatorEverywhere : public Indicator
{   

public:

   /// Constructor
   IndicatorEverywhere
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   ) : Indicator(_avgr, _fes, _fes_const, _offsets, _d) {};

   /// Destructor
   ~IndicatorEverywhere() {};

   /// Find troubled cells
   virtual double checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) 
   { 
      setValue(iCell, 0.0); return 0.0;
   };
};

#endif // INDICATOR_EVERYWHERE_H

