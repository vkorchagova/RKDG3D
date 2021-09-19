#ifndef STENCIL_H
#define STENCIL_H

#include "mfem.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

/// Proc rank 
extern int myRank;


class Stencil
{

public:

   /// Stencil for current cell (cell numbers)
   Array<int> cell_num;

   /// Max stencil size
   const int max_size = 10;

   /// Internal face numbers in stencil
   Array<int> internal_face_numbers;

   /// Shared face numbers in stencil
   Array<int> shared_face_numbers;

public:

   Stencil()
   {
      cell_num.Reserve(max_size);
      internal_face_numbers.Reserve(max_size);
      shared_face_numbers.Reserve(max_size);
   };

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