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
   LimiterFinDiff(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace *_fes, const Array<int>& _offsets, int _d) 
      : Limiter(_ind,_avgr,_fes,_offsets,_d) {};
   
   /// Destructor
   ~LimiterFinDiff() {};

   /// Limit solution
   virtual void limit(const int iCell, const Vector& el_ind, DenseMatrix& elfun1_mat) override;  
};

#endif // LIMITER_FINDIFF_H