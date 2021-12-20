#include "limiter_findiff.hpp"

void LimiterFinDiff::limit(const int iCell, const Vector& el_ind, DenseMatrix& elfun1_mat) 
{
   bool totalInd = false;

   for (int iEq = 0; iEq < num_equation; ++iEq)
      totalInd = el_ind[iEq] < 0.999 ? true : false;

   // get FE
   fe = fes->GetFE(iCell);

   const int nDofs = fe->GetDof();

   averager.readElementAverageByNumber(iCell, el_uMean);

   // if (iCell == 2703 || iCell == 1191) 
   // {
   //    el_uMean.Print(cout << iCell << ": el_uMean = ");
   // }

   // replace solution to mean values
   for (int iEq = 0; iEq < num_equation; ++iEq)
      for (int iDof = 0; iDof < nDofs; ++iDof)
      {
         // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval before lim = " << elfun1_mat(iDof, iEq) << endl;}
         elfun1_mat(iDof, iEq) = totalInd ? el_uMean(iEq) : elfun1_mat(iDof, iEq);
         // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval after lim = " << elfun1_mat(iDof, iEq) << endl;}
      }

      // // std::cout << "=======\n el_x inside limiter = ";
      // //    el_x.Print(std::cout); 
};   
