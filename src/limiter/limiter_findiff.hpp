#ifndef LIMITER_FINDIFF_H
#define LIMITER_FINDIFF_H

#include "limiter.hpp"

/// 
/// Simple limiting algorithm with slope zerofication
///
class LimiterFinDiff : public Limiter
{
public:

   /// Constructor 
   LimiterFinDiff(ParFiniteElementSpace *_fes, const Array<int>& _offsets, int _d) 
      : Limiter(_fes,_offsets,_d) {};
   
   /// Destructor
   ~LimiterFinDiff() {};

   /// Limit solution
   virtual void limit(Vector &x) override;  
};

#endif // LIMITER_FINDIFF_H