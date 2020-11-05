#include "limiter_findiff.hpp"

void LimiterFinDiff::limit(Vector &x) 
{ 
   // std:: cout << "FinDiff\n";

   ParMesh* mesh = fes->GetParMesh();

   /// Mean value of solution in one cell
   Vector uMean(num_equation);  

   const FiniteElement *fe1;

   // loop through finite elements
   // cout << "par mesh NE = " << mesh->GetNE() << endl;
   // cout << "x size = " << x.Size() << endl;
   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
      // std::cout << iCell << std::endl;

      // get FE
      fe1 = fes->GetFE(iCell);

      // get FE dofs indices
      fes->GetElementVDofs(iCell, vdofs);
      const int nDofs = fe1->GetDof();

      //get FE transformation
      el_trans = mesh->GetElementTransformation(iCell);
      //el_x.SetSize(nDofs);

      // read sol values for defined vdofs
      x.GetSubVector(vdofs, el_x);

      // make matrix from values data
      DenseMatrix elfun1_mat(el_x.GetData(), nDofs, num_equation);
      // std::cout << "elfun1_mat = ";
      // elfun1_mat.Print(std::cout);

      // get center of finite element
      mesh->GetElementCenter(iCell, el_center_phys);
      //if (iCell == 0) {std::cout << "el_center_phys = "; el_center_phys.Print(std::cout);}

      // move it to the element space
      el_trans->TransformBack(el_center_phys, el_center_ref);
      // std::cout << "el_center_ref = " << el_center_ref.index << " " << el_center_ref.x << " " << el_center_ref.y << " " << el_center_ref.z << std::endl;

      // get values of shape functions in cell center
      shape.SetSize(nDofs);
      fe1->CalcShape(el_center_ref, shape);
      // std::cout << "shape = ";
      // shape.Print(std::cout);

      // restore the mean value of solution (it is solution in center for linear function)
      elfun1_mat.MultTranspose(shape, uMean);
      // std::cout << "uMean = ";
      // uMean.Print(std::cout);

      // replace solution to mean values
      for (int iEq = 0; iEq < num_equation; ++iEq)
         for (int iDof = 0; iDof < nDofs; ++iDof)
         {
            // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval before lim = " << elfun1_mat(iDof, iEq) << endl;}
            elfun1_mat(iDof, iEq) = uMean(iEq);
            // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval after lim = " << elfun1_mat(iDof, iEq) << endl;}
         }

      // // std::cout << "=======\n el_x inside limiter = ";
      // //    el_x.Print(std::cout);

      x.SetSubVector(vdofs, el_x);

   
   }  
   // std::cout << "=======\n X inside limiter = ";
   // x.Print(std::cout); 
};   
