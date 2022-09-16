#ifndef INDICATOR_MICHALAK_H
#define INDICATOR_MICHALAK_H

#include "indicator_bj.hpp"

/// 
/// Barth-Jespersen indicator based on limiting algorithm
///
class IndicatorMichalak : public IndicatorBJ
{

private:

   /// Compute correction function value
   virtual void computeCorrectionFunction(const double& y, double& yMin);
  

public:
   
   /// Constructor
   IndicatorMichalak
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   );

   /// Destructor
   ~IndicatorMichalak() {};

};

#endif // INDICATOR_MICHALAK_H