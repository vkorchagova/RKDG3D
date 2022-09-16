#ifndef STENCIL_H
#define STENCIL_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

///
/// Class for cell stencil needed for limiters
///
class Stencil
{

public:

   /// Stencil for current cell (cell numbers)
   Array<int> cell_num;

   /// Max stencil size (just for reserve)
   const int max_size = 10;

   /// Internal face numbers in stencil
   Array<int> internal_face_numbers;

   /// Shared face numbers in stencil
   Array<int> shared_face_numbers;

public:

   /// Constructor
   Stencil()
   {
      cell_num.Reserve(max_size);
      internal_face_numbers.Reserve(max_size);
      shared_face_numbers.Reserve(max_size);
   };

   /// Destructor
   ~Stencil() {};

   // Clean stencil values
   void clean()
   {
      internal_face_numbers.DeleteAll();
      shared_face_numbers.DeleteAll();
      cell_num.DeleteAll();
   };

   /// Print values
   void Print(ParMesh& mesh, int iCell, int rank = 0)
   {
      cout << "======" << endl;
      cout << "Stencil for cell #" << iCell << ", proc #" << rank << endl;
      cout << "------" << endl;
      cell_num.Print(cout << "cell num = ");
      internal_face_numbers.Print(cout << "internal_face_numbers num = ");
      shared_face_numbers.Print(cout << "shared_face_numbers num = ");
      cout << "Cell centres:" << endl;
      for (int iCellNum : cell_num)
      {
         // int iCellNumOk = iCellNum;
         // if (iCellNum > mesh->GetNE())
         //    iCellNumOk(avgs->FaceNbrData())[iEq  + (iCell - mesh->GetNE()) * num_equation];
         if (iCellNum > mesh.GetNE())
         {
            FaceElementTransformations* face_el_trans = mesh.GetSharedFaceTransformations(shared_face_numbers[0]);
            cout << face_el_trans->Elem1No << ' ' << face_el_trans->Elem2No << endl;
         }
         else
         {
            Vector cellCenter(3);
            mesh.GetElementCenter(iCellNum,cellCenter);
            cellCenter.Print(cout << "\t centre of cell #" << iCellNum);
         }
      }
      cout << "======" << endl;
   }
};

#endif // STENCIL_H