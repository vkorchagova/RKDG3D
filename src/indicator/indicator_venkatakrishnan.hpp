#ifndef INDICATOR_VENKATAKRISHNAN_H
#define INDICATOR_VENKATAKRISHNAN_H

#include "indicator_bj.hpp"

/// 
/// Barth-Jespersen indicator based on limiting algorithm
///
class IndicatorVenkatakrishnan : public IndicatorBJ
{
private:

   /// Compute correction function value
   virtual void computeCorrectionFunction(const double& y, double& yMin);
  

public:
   
   /// Constructor
   IndicatorVenkatakrishnan
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   );

   /// Destructor
   ~IndicatorVenkatakrishnan() {};

};

#endif // INDICATOR_VENKATAKRISHNAN_H