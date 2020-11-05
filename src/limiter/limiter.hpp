#ifndef LIMITER_H
#define LIMITER_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern const int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for limiting algorithm
/// for the DG slopes
///
class Limiter
{
protected:

   /// Mesh
   ParMesh* mesh;

   /// FE finite element space (connection with mesh and dofs)
   ParFiniteElementSpace *fes;

   /// Array of element DOF indices
   Array<int> vdofs;

   /// Solution in element dofs 
   Vector el_x;

   /// Center of FE in physical space
   Vector el_center_phys;

   /// Center of FE in reference space
   IntegrationPoint el_center_ref;

   /// Physical-to-reference and vice versa transformation rule for FE
   ElementTransformation* el_trans;
   
   /// Values of shape functions in defined point
   Vector shape;

   /// Space dimension
   int dim;

   /// parGridFunction to enable special pgrf functions for x vector
   mutable ParGridFunction parGridX;

   /// Offsets to deal with variables component-by-component
   const Array<int>& offsets;

   /// Compute average of extrapolated function values from neighbours on troubled cell
   void GetStencilAverages(
      const GridFunction& x,
      const Array<int>& stencil_num,
      GridFunction& avgs
   ) const;

   /// Special function to compute shape fun values on troubled cell via neighbour cell
   void AssembleShiftedElementMatrix(
      const FiniteElement &trial_fe, 
      const FiniteElement &troubled_fe,
      const FiniteElement &test_fe,
      ElementTransformation &Trans, 
      DenseMatrix &elmat
   ) const;

public:

   /// Constructor
   Limiter(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d) 
      : fes(_fes), offsets(_offsets), dim(_d) {parGridX.MakeRef(_fes, NULL);};

   /// Destructor
   virtual ~Limiter() {};

   /// Limit solution
   virtual void limit(Vector &x) = 0;
   
};

#endif // LIMITER_H

