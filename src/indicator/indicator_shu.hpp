#ifndef INDICATOR_SHU_H
#define INDICATOR_SHU_H

#include "indicator.hpp"

/// 
/// Shu indicator
///
class IndicatorShu : public Indicator
{
private:

   /// Indication constant
   const double Ck = 0.03;

   /// Unzero-denom
   double eps = 1e-6;

   /// max | u_avg_j|, j \in stensil
   Vector maxFabsPj;

   /// sum | u_avg_0 - u_avg_extrap_j|, j \in stensil_neibs
   Vector sumFabsDiffExtrap;
  

public:
   
   /// Constructor
   IndicatorShu(Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata);
   
   /// Destructor
   ~IndicatorShu()
   {
      
   };

   /// Limit solution
   virtual void checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) override;
};

#endif // INDICATOR_SHU_H
