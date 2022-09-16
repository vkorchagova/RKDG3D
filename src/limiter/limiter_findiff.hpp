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
   LimiterFinDiff(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace *_fes, const Array<int>& _offsets, bool _linearize = 0, bool _haveLastHope = 1, int _fdGroupAttribute = -1, int _d = 3) 
      : Limiter(_ind,_avgr,_fes,_offsets,_linearize,_haveLastHope,_fdGroupAttribute,_d) {};
   
   /// Destructor
   ~LimiterFinDiff() {};

   /// Limit solution
   virtual void limit(const int iCell, const double ind_value, const double nDofs, DenseMatrix& elfun1_mat) override;  
};

#endif // LIMITER_FINDIFF_H