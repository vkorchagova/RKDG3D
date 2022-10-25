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
   double mI;

   /// Maximal averaged value in stencil
   double MI;

   /// Average values on one element
   double el_uMean;

   /// Indicator value
   double yMin;

   /// Index of solution component to be chosen as "root" for computing correction function
   const int iSolRoot = DEFAULT_BJ_ROOT_SOL_COMP_INDEX;

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
