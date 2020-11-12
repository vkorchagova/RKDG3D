#ifndef INDICATOR_SHU_H
#define INDICATOR_SHU_H

#include "indicator.hpp"

/// 
/// Barth-Jespersen indicator based on limiting algorithm
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

   /// Values u_avg_extrap_j
   ParGridFunction* avgs_extrap;
   BlockVector* u_block_avg_extrap;
   Array<int> offsets_avg_extrap;
   ParFiniteElementSpace *fes_avg_extrap;
   DG_FECollection* fec_avg_extrap;
   ParFiniteElementSpace* fes_avg_extrap_component;
  

public:
   
   /// Constructor
   IndicatorShu(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata);
   
   /// Destructor
   /// Destructor
   ~IndicatorShu()
   {
      // delete stencil;
      delete avgs_extrap; 
      delete u_block_avg_extrap; 
      delete fes_avg_extrap_component;
      delete fes_avg_extrap; 
      delete fec_avg_extrap;
   };

   /// Limit solution
   virtual void checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const ParGridFunction* uMean, 
      const DenseMatrix& elfun1_mat,
      ParGridFunction &x) override;
};

#endif // INDICATOR_SHU_H
