#ifndef AVERAGER_H
#define AVERAGER_H

#include "mfem.hpp"
#include "stencil.hpp"


using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Class to save average values in cells
/// ... and average values of extrapolated solution in cell neighbours
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


   /// Offsets to deal with variables component-by-component
   const Array<int>& offsets;

   /// Actual solution data
   Vector* x;

   /// Actual solution data as a ParGridFunction
   ParGridFunction* parGridX;


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


   /// Finite element collection for average values of full solution
   DG_FECollection* fec_avg;

   /// Finite element space for average values of full solution
   ParFiniteElementSpace *fes_avg;

   /// Finite element space for average values of one component of solution
   ParFiniteElementSpace* fes_avg_component;
   
   /// Block vector to store average values
   BlockVector* u_block_avg;

   /// Offsets for looping through average values
   Array<int> offsets_avg;

   /// ParGridFunction to deal with average values
   ParGridFunction* avgs;


   /// Finite element collection for average values of full extrapolated solution to neighbour
   DG_FECollection* fec_avg_extrap;

   /// Finite element space for average values of full extrapolated solution to neighbour
   ParFiniteElementSpace *fes_avg_extrap;

   /// Finite element space for average values of one component of extrapolated solution to neighbour
   ParFiniteElementSpace* fes_avg_extrap_component;

   /// BlockVector to store with average values of extrapolated solution to neighbour
   BlockVector* u_block_avg_extrap;

   /// Offsets for looping through average values of extrapolated solution to neighbour
   Array<int> offsets_avg_extrap;

   /// ParGridFunction to deal with average values of extrapolated solution to neighbour
   ParGridFunction* avgs_extrap;


   /// Special function to compute shape function values on troubled cell via neighbour cell
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
   ~Averager();

   /// Update actual solution values
   void update( Vector* _x, ParGridFunction* _parGridX) { x = _x; parGridX = _parGridX; };

   /// Compute mean values in all cells in domain via honest GridFunction
   void computeMeanValues();

   /// Read element average: special function due to different mechanisms for internal and shared values
   void readElementAverageByNumber(
      const int iCell, 
      Vector& el_uMean
   );
   void readElementAverageComponentByNumber(
      const int iCell, 
      const int iEq, 
      double value
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

   /// Compute one component of average of extrapolated function values from neighbours on troubled cell
   void computeStencilExtrapAverages(
      const ParGridFunction& x,
      const Stencil* stencil,
      ParGridFunction& avgs
   );  


   /// Update FESpaces when mesh is changed (for AMR)
   void updateSpaces();

   /// Update all grid functions when mesh is changed (for AMR)
   void updateSolutions();

   /// Send finish of update when mesh is changed (for AMR)
   void updateFinished();
};



#endif // AVERAGER_H

