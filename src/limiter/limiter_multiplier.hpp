#ifndef LIMITER_MULTIPLIER_H
#define LIMITER_MULTIPLIER_H

#include "limiter.hpp"

/// 
/// Limitation by simple multiplication of slopes to predefined coefficients
///
class LimiterMultiplier : public Limiter
{
   /// Values of shape functions in defined point
   Vector el_shape;

   /// Additional linearization (zerofication of XY in rectangular and hexahedral cells)
   void linearize(const int iCell, const Vector& el_uMean, DenseMatrix &elfun1_mat);

public:
   
   /// Constructor
   LimiterMultiplier(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace *_fes, const Array<int>& _offsets, bool _linearize = 0, bool _haveLastHope = 1, int _fdGroupAttribute = -1, int _d = 3);
   
   /// Destructor
   ~LimiterMultiplier() {};

   /// Limit solution
   virtual void limit(const int iCell, const double ind_value, DenseMatrix& elfun1_mat) override;
};

#endif // LIMITER_MULTIPLIER_H
