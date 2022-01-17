#ifndef LIMITER_NONE_H
#define LIMITER_NONE_H

#include "limiter.hpp"

/// 
/// Just return solution without any limitation
///
class LimiterNone : public Limiter
{

public:

   /// Constructor 
   LimiterNone(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace *_fes, const Array<int>& _offsets, bool _linearize = 0, bool _haveLastHope = 1, int _fdGroupAttribute = -1, int _d = 3) 
      : Limiter(_ind,_avgr,_fes,_offsets,_linearize,_haveLastHope, _fdGroupAttribute, _d) {};
   
   /// Destructor
   ~LimiterNone() {};

   /// Limit solution
   virtual void limit(const int iCell, const double ind_value, DenseMatrix& elfun1_mat) override {};  
};

#endif // LIMITER_NONE_H