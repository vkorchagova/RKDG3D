#include "indicator_bj.hpp"
 
// bool pleaseWriteMe = false;

IndicatorBJ::IndicatorBJ
(
   Averager& _avgr, 
   ParFiniteElementSpace* _fes,
   ParFiniteElementSpace* _fes_const, 
   const Array<int>& _offsets, 
   int _d
) : Indicator(_avgr, _fes, _fes_const, _offsets, _d) 
{
   mI.SetSize(num_equation);
   MI.SetSize(num_equation);
   el_uMean.SetSize(num_equation);
   yMin.SetSize(num_equation);
};

void IndicatorBJ::computeCorrectionFunction(const double& y, double& yMin)
{
   yMin = y < yMin ? y : yMin;
   yMin = yMin > 1.0 ? 1.0 : yMin;
}

void IndicatorBJ::updateYmin(
   const IntegrationRule& ir, 
   IntegrationPointTransformation* curTrans, 
   const DenseMatrix& elfun1_mat, 
   const int iCell) 
{
   IntegrationPoint eip1;
   Vector funval1(num_equation);
   Vector diff(num_equation);
   Vector y(num_equation);

   for (int iPoint = 0; iPoint < ir.GetNPoints(); iPoint++)
   {
      const IntegrationPoint &ip = ir.IntPoint(iPoint);

      if (curTrans)
         curTrans->Transform(ip, eip1);
      else
         eip1 = ip;

      // Calculate basis functions on elements at the face
      el_shape.SetSize(fe->GetDof());
      fe->CalcShape(eip1, el_shape);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(el_shape, funval1);


      averager.readElementAverageByNumber(iCell, el_uMean);
      subtract(funval1, el_uMean, diff);

      // Compute different fractions
      for (int i = 0; i < num_equation; ++i)
      // for (int i = 0; i < 1; ++i)
      {
         if (diff[i] > DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean[i])))
               y[i] = (MI[i] - el_uMean[i]) / diff[i];
          else if (diff[i] < - DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean[i])))
               y[i] = (mI[i] - el_uMean[i]) / diff[i];
          else
               y[i] = 1.0;

         // if (iCell == 1512 || iCell == 0)
         // {
         //   std::cout  << std::setprecision(18)
         //         << " * iCell = " << iCell
         //         << " -- y[i] = " << y[i] 
         //         << "\n -- MI[i] = " << MI[i]
         //         << "\n -- mI[i] = " << mI[i]
         //         << "\n -- funval1[i] = " << funval1[i]
         //         << "\n -- el_uMean[i] = " << el_uMean[i]
         //         << "\n -- diff[i] = " << diff[i]
         //          << std::setprecision(6)
         //         << std::endl;
         // }
         computeCorrectionFunction(y[i],yMin[i]);
      
          
          if (yMin[i] < 0)
               std::cout << " yMin = " << yMin[i] << " < 0 --- so strange! " << MI[i] << ' ' << mI[i] << ' ' << funval1[i] << ' ' << el_uMean[i] << ' ' << diff[i] << ' ' << mI[i] - el_uMean[i]  << ' ' << MI[i] - el_uMean[i] << "; iCell = " << iCell << "; i = " << i << std::endl;
      }
   } // for iPoint
}

double IndicatorBJ::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil, 
   const DenseMatrix& elfun1_mat)
{
   // loop through finite elements to compute indicator values

   mI = DEFAULT_LARGE_NUMBER;
   MI = -DEFAULT_LARGE_NUMBER;
   yMin = DEFAULT_LARGE_NUMBER;

   // get FE
   fe = fes->GetFE(iCell);
   fes->GetElementVDofs(iCell, el_vdofs);
   const int nDofs = fe->GetDof();

   // if (iCell == 1512 || iCell == 0)
   // { 
   //   std::cout << "iCell = " << iCell << std::endl;
   //   std::cout << "elfun1_mat = ";
   //   elfun1_mat.Print(std::cout << std::setprecision(18) );
   // }

   // std::cout << "stencil: ";
   // stencil_num.Print(std::cout);

   // compute min|max values of solution for current stencil
   
   
   for (int k : stencil->cell_num)
   {
      // if (iCell == 3160) 
      // {
      //   el_uMean.Print(std::cout << "el_uMean for cell #" << k << ":");
      // }

      averager.readElementAverageByNumber(k, el_uMean);
      for (int iSol = 0; iSol < num_equation; ++iSol)
      {
         
         mI[iSol] = el_uMean[iSol] < mI[iSol] ? el_uMean[iSol] : mI[iSol];
         MI[iSol] = el_uMean[iSol] > MI[iSol] ? el_uMean[iSol] : MI[iSol];
      }
   }

   // if (iCell == 1789 && myRank == 17) 
   // {
   //   std::cout << "stencil = ";
   //   stencil_num.Print(std::cout);
   //   std::cout << "internal faces = ";
   //   internal_face_numbers.Print(std::cout);
   // }


   // run through faces again to compute values in gauss points
   for (int iFace : stencil->internal_face_numbers)
   {
      // Integration order calculation from DGTraceIntegrator
      face_el_trans = mesh->GetFaceElementTransformations(iFace);
      IntegrationPointTransformation curTrans = face_el_trans->Elem1No == iCell ? face_el_trans->Loc1 : face_el_trans->Loc2;

      int intorder = 2;//(std::min(face_el_trans->Elem1->OrderW(), face_el_trans->Elem2->OrderW()) +
                      //2*std::max(fes->GetFE(face_el_trans->Elem1No)->GetOrder(), fes->GetFE(face_el_trans->Elem2No)->GetOrder()));

      const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);
      // if (iCell == 3160) 
      // {
      //   std::cout << "tr->FaceGeom = " << face_el_trans->FaceGeom << ", intorder = " << intorder << std::endl;
      //   std::cout << "npoints = " << ir->GetNPoints() << std::endl;
      //   std::cout << "numFace = " << iFace << ";\n";
      // }

      updateYmin(*ir, &curTrans, elfun1_mat, iCell);

   } // for iFace

   // run through faces again to compute values in gauss points
   for (int iFace : stencil->shared_face_numbers)
   {
      // Integration order calculation from DGTraceIntegrator
      face_el_trans = mesh->GetSharedFaceTransformations(iFace);
      IntegrationPointTransformation curTrans = face_el_trans->Loc1;

      int intorder;
         intorder = face_el_trans->Elem1->OrderW() + 2*fes->GetFE(face_el_trans->Elem1No)->GetOrder();

      const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);

      // if (iCell == 1789 && myRank == 17) 
      // {
      //   std::cout << "tr->FaceGeom = " << face_el_trans->FaceGeom << ", intorder = " << intorder << std::endl;
      //   std::cout << "npoints = " << ir->GetNPoints() << std::endl;
      //   std::cout << "numFace = " << iFace << ";\n";
      // }
      updateYmin(*ir, &curTrans, elfun1_mat, iCell);
   } // for iFace

   // if (iCell == 1512 || iCell == 0) 
   // {
   //   yMin.Print(std::cout << "yMin = ");
   // }

   // double iVal = 1e+6;
   // for (double a : yMin)
   //   if (iVal > a)
   //       iVal = a;



   // set indicator values to the external field 
   double iVal = yMin[0];//int(yMin[0]*1e+8 + 0.5)/1e+8;//yMin[0];// > 0.999999 ? 1.0 : yMin[0];
 
   setValue(iCell, iVal);

   return iVal;
};
