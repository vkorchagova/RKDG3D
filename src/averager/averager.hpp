#ifndef AVERAGER_H
#define AVERAGER_H

#include "mfem.hpp"
#include "stencil.hpp"

using namespace std;
using namespace mfem;

extern const int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for limiting algorithm
/// for the DG slopes
///
class Averager
{
protected:

   /// Mesh
   ParMesh* mesh;

   /// FE finite element space (connection with mesh and dofs)
   ParFiniteElementSpace *fes; 

   /// Space dimension
   int dim;


   /// Actual solution data
   Vector* x;
   ParGridFunction* parGridX;

   /// Offsets to deal with variables component-by-component
   const Array<int>& offsets;





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

   

   


   /// Finite element space for average values
   DG_FECollection* fec_avg;
   ParFiniteElementSpace *fes_avg;
   

   /// Mean values of solution in all cells
   /// for high order vals
   ParGridFunction* avgs;
   BlockVector* u_block_avg;

   /// Offsets for looping through average values
   Array<int> offsets_avg;


   /// Values u_avg_extrap_j
   ParGridFunction* avgs_extrap;
   BlockVector* u_block_avg_extrap;
   Array<int> offsets_avg_extrap;
   ParFiniteElementSpace *fes_avg_extrap;
   DG_FECollection* fec_avg_extrap;
   ParFiniteElementSpace* fes_avg_extrap_component;

   /// Special function to compute shape fun values on troubled cell via neighbour cell
   void assembleShiftedElementMatrix(
      const FiniteElement &trial_fe, 
      const FiniteElement &troubled_fe,
      const FiniteElement &test_fe,
      ElementTransformation &Trans, 
      ElementTransformation &TransTroubled,
      DenseMatrix &elmat
   ); 

public:

   /// Constructor
   Averager(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d); 

   /// Destructor
   virtual ~Averager()
   {
      delete avgs_extrap; 
      delete u_block_avg_extrap; 
      delete fes_avg_extrap_component;
      delete fes_avg_extrap; 
      delete fec_avg_extrap;

      delete avgs; 
      delete u_block_avg; 
      delete fes_avg; 
      delete fec_avg;
   };

   /// Update actual solution values
   void update( Vector* _x, ParGridFunction* _parGridX) { x = _x; parGridX = _parGridX; };

   /// Compute mean values in all cells in domain via honest GridFunction
   void computeMeanValues();

   /// Read element average: special function due to different mechanisms for internal and shared values
   void readElementAverageByNumber(
      const int iCell, 
      Vector& el_uMean
   );

   /// Read element extrap average: special function due to different mechanisms for internal and shared values
   void readElementExtrapAverageByNumber(
      const int iCell, 
      Vector& el_uMean
   );

   /// Compute average of extrapolated function values from neighbours on troubled cell
   void computeStencilExtrapAveragesVector(
      const Stencil* stencil
   );

   void computeStencilExtrapAverages(
      const ParGridFunction& x,
      const Stencil* stencil,
      ParGridFunction& avgs
   );  
};



#endif // AVERAGER_H

