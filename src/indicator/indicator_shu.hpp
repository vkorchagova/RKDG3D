#ifndef INDICATOR_SHU_H
#define INDICATOR_SHU_H

#include "indicator.hpp"

/// 
/// Shu indicator
///
class IndicatorShu : public Indicator
{
private:

   /// max | u_avg_j|, j \in stensil
   Vector maxFabsPj;

   /// sum | u_avg_0 - u_avg_extrap_j|, j \in stensil_neibs
   Vector sumFabsDiffExtrap;
  

public:
   
   /// Constructor
   IndicatorShu
   (
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   );

   /// Destructor
   ~IndicatorShu() {};

   /// Limit solution
   virtual double checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) override;
};

#endif // INDICATOR_SHU_H
