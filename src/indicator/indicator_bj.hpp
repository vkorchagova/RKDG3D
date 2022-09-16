#ifndef INDICATOR_BJ_H
#define INDICATOR_BJ_H

#include "indicator.hpp"

/// 
/// Barth-Jespersen indicator based on limiting algorithm
///
class IndicatorBJ : public Indicator
{
private:

   /// Minimal averaged value in stencil
   Vector mI;

   /// Maximal averaged value in stencil
   Vector MI;

   /// Average values on one element
   Vector el_uMean;

   /// Indicator value
   Vector yMin;

   /// Update limiter coefficient
   void updateYmin(
      const IntegrationRule& ir, 
      IntegrationPointTransformation* curTrans, 
      const DenseMatrix& elfun1_mat,
      const int iCell 
   );

   /// Compute correction function value
   virtual void computeCorrectionFunction(const double& y, double& yMin);
  

public:
   
   /// Constructor
   IndicatorBJ
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   );

   /// Destructor
   ~IndicatorBJ() {};

   /// Limit solution
   virtual double checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) override;
};

#endif // INDICATOR_BJ_H
