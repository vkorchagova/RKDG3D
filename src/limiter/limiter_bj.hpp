#ifndef LIMITER_BJ_H
#define LIMITER_BJ_H

#include "limiter.hpp"

/// 
/// Barth-Jespersen limiting algorithm
///
class LimiterBJ : public Limiter
{
private:

   /// Use primitive variables
   const bool primitive = false;

   /// Mean values of solution in all cells
   DenseMatrix uMean;

   /// Stencil for current cell
   Array<int> stencil_num;

   /// Pointer to finite element
   const FiniteElement *fe;

   /// Vector to store limited values
   Vector xNew;

   const int stencil_max_size = 10;

   /// Minimal averaged value in stencil
   Vector mI;

   /// Maximal averaged value in stencil
   Vector MI;

   /// Transformation rule on face
   FaceElementTransformations* tr;



   //....... the foolowing vars only for future

   /// Stencils for each cell
   Array<Array<int>> stencil;

   /// Offsets for looping through average values
   Array<int> offsets_avg;

   /// Finite element space for average values
   ParFiniteElementSpace *fes_avg;

   /// Mean values of solution in all cells
   /// for high order vals
   ParGridFunction* avgs;
   BlockVector* u_block_avg;
   DG_FECollection* fec_avg;


   /// Initialisation: find stencils for all cells
   void Init(const ParMesh* mesh);

   /// Compute element average
   void linearize(const int iCell, const DenseMatrix & uMean, DenseMatrix &elfun1_mat);

   /// Update limiter coefficient
   void updateYMin(
      const IntegrationRule& ir, 
      IntegrationPointTransformation* curTrans, 
      const DenseMatrix& elfun1_mat, 
      const int iCell, 
      Vector &yMin 
   );
  

public:
   
   /// Constructor
   LimiterBJ(ParFiniteElementSpace *_fes, const Array<int>& _offsets, int _d);
   
   /// Destructor
   ~LimiterBJ()
   {
      delete avgs; 
      delete u_block_avg; 
      delete fes_avg; 
      delete fec_avg;
   };

   

   /// Limit solution
   virtual void limit(Vector &x) override;
};

#endif // LIMITER_BJ_H
