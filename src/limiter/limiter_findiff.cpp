#include "limiter_findiff.hpp"

void LimiterFinDiff::limit(const int iCell, const double ind_value, const double nDofs, DenseMatrix& elfun1_mat) 
{
   averager.readElementAverageByNumber(iCell, el_uMean);

   // if (iCell == 2703 || iCell == 1191) 
   // {
   //   el_uMean.Print(std::cout << iCell << ": el_uMean = ");
   // }

   // replace solution to mean values
   for (int iEq = 0; iEq < num_equation; ++iEq)
      for (int iDof = 0; iDof < nDofs; ++iDof)
      {
         // if (myRank == 0 && iCell == 0 && iEq == 0) {std::cout << "   funval before lim = " << elfun1_mat(iDof, iEq) << std::endl;}
         elfun1_mat(iDof, iEq) = ind_value < 0.999999 ? el_uMean(iEq) : elfun1_mat(iDof, iEq);
         // if (myRank == 0 && iCell == 0 && iEq == 0) {std::cout << "   funval after lim = " << elfun1_mat(iDof, iEq) << std::endl;}
      }

      // // std::cout << "=======\n el_x inside limiter = ";
      // //   el_x.Print(std::cout); 
};   
