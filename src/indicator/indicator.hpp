#ifndef INDICATOR_H
#define INDICATOR_H

#include "mfem.hpp"
#include "averager.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for search of cells where solution should be limited
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

   /// FE finite element space for constant-piecewise elements
   ParFiniteElementSpace *fes_const;

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

   /// Values of indicator field associated with ParaView external writer
   ParGridFunction* values;
   
public:
   
   ParGridFunction* getValues() { return values; };

   /// Constructor
   Indicator(
      Averager& _avgr, 
      ParFiniteElementSpace* _fes,
      ParFiniteElementSpace* _fes_const, 
      const Array<int>& _offsets, 
      int _d
   );

   /// Destructor
   virtual ~Indicator();

   /// Check if defined cell is troubled and set indicator value
   virtual void checkDiscontinuity(
      const int iCell, 
      const Stencil* stencil, 
      const DenseMatrix& elfun1_mat
   ) = 0;

   /// Set indicator value for defined cell
   void setValue(int iCell, double val);

   /// Read indicator value for defined cell
   const double getValue(int iCell);


};

#endif // INDICATOR_H

