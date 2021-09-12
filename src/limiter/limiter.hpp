#ifndef LIMITER_H
#define LIMITER_H

#include "mfem.hpp"
#include "indicator.hpp"
#include "averager.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

/// Proc rank 
extern int myRank;

class Stencil;
class Indicator;

/// 
/// Abstract class for limiting algorithm
/// for the DG slopes
///
class Limiter
{
protected:

   /// Averager
   Averager& averager;

   /// Mesh
   ParMesh* mesh;

   /// Indicator
   Indicator& indicator;

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
   
   

   

   /// parGridFunction to enable special pgrf functions for x vector
   mutable ParGridFunction parGridX;

   /// Offsets to deal with variables component-by-component
   const Array<int>& offsets;

   /// Stencil
   Stencil* stencil;   

   /// Average values on one element
   Vector el_uMean;

   /// Loc-to-global mapping for shared faces
   std::map<int,int> lf2sf;

   /// Read element average due to different mechanisms for internal and shared values
   friend void readElementAverageByNumber(const int iCell, Vector& el_uMean);

   /// Compute mean values via honest GridFunction
   void computeMeanValues();

   /// Compute stencil
   void getStencil(const int iCell);
   void cleanStencil();



public:

   /// Constructor
   Limiter(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d); 

   /// Destructor
   virtual ~Limiter()
   {
      delete stencil;
   };

   /// Update solution (wrapper for checkDiscontinuities + limit)
   void update(Vector &x);

   /// Limit solution
   virtual void limit(const int iCell, const Vector& el_ind, DenseMatrix& elfun1_mat) = 0;
   
};

#endif // LIMITER_H

