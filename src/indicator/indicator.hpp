#ifndef INDICATOR_H
#define INDICATOR_H

#include "mfem.hpp"
#include "averager.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for checking discontinuities
/// for the DG slopes
///
class Indicator
{

protected:

   /// Averager
   Averager& averager;

   /// Mesh
   ParMesh* mesh;

   /// FE finite element space (connection with mesh and dofs)
   ParFiniteElementSpace *fes;

   /// Space dimension
   int dim;

   /// Vector to store limited values
   Vector xNew;


   /// Pointer to finite element
   const FiniteElement *fe;

   /// DOF indices for element
   Array<int> el_vdofs;

   /// Solution in element dofs 
   Vector el_x;

   /// Physical-to-reference and vice versa transformation rule for FE
   ElementTransformation* el_trans;

   /// Transformation rule from face space to element space
   FaceElementTransformations* face_el_trans;
   
   /// Values of shape functions in defined point
   Vector el_shape;

   /// parGridFunction to enable special pgrf functions for x vector
   mutable ParGridFunction parGridX;

   /// Offsets to deal with variables component-by-component
   const Array<int>& offsets;
   

public:

   /// Values of indicator field associated with ParaView external writer
   BlockVector& values;

   /// Keep minimal indicator value 
   Vector minValues;

   /// Constructor
   Indicator(Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata); 

   /// Destructor
   virtual ~Indicator() {};

   /// Find troubled cells
   virtual void checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) = 0;
};

#endif // INDICATOR_H

