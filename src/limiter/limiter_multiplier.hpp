#ifndef LIMITER_MULTIPLIER_H
#define LIMITER_MULTIPLIER_H

#include "limiter.hpp"

/// 
/// Barth-Jespersen limiting algorithm
///
class LimiterMultiplier : public Limiter
{
   /// Additional linearization (zerofication of XY in rectangular and hexahedral cells)
   void linearize(const int iCell, const Vector& el_uMean, DenseMatrix &elfun1_mat);

public:
   
   /// Constructor
   LimiterMultiplier(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace *_fes, const Array<int>& _offsets, int _d);
   
   /// Destructor
   ~LimiterMultiplier() {};

   /// Limit solution
   virtual void limit(const int iCell, const Vector& el_ind, DenseMatrix& elfun1_mat) override;
};

#endif // LIMITER_MULTIPLIER_H
