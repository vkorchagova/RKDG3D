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
};

#endif // STENCIL_H