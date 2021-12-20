#include "indicator_bj.hpp"

// bool pleaseWriteMe = false;

IndicatorBJ::IndicatorBJ(Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata)
    : Indicator(_avgr, _fes, _offsets, _d, _idata) 
{
   mI.SetSize(num_equation);
   MI.SetSize(num_equation);
   el_uMean.SetSize(num_equation);
   yMin.SetSize(num_equation);
};


void IndicatorBJ::updateYMin(
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
      {
         if (diff[i] > 1e-3)// * max(1.0, fabs(el_uMean[i])))
              y[i] = (MI[i] - el_uMean[i]) / diff[i];
          else if (diff[i] < - 1e-3)// * max(1.0, fabs(el_uMean[i])))
              y[i] = (mI[i] - el_uMean[i]) / diff[i];
          else
              y[i] = 1.0;

         // if (iCell == 3160 && i == 0)
         // {
         //    cout << " -- y[i] = " << y[i] 
         //         << "\n -- MI[i] = " << MI[i]
         //         << "\n -- mI[i] = " << mI[i]
         //         << "\n -- el_uMean[i] = " << el_uMean[i]
         //         << "\n -- diff[i] = " << diff[i]
         //         << endl;
         // }
         
         // // classical BJ
         // yMin[i] = y[i] < yMin[i] ? y[i] : yMin[i];
         // yMin[i] = yMin[i] > 1.0 ? 1.0 : yMin[i];

         // // venkatakrishnan
         // double yCur = (y[i]*y[i] + 2.0 * y[i]) / (y[i] * y[i] + y[i] + 2.0);
         // yMin[i] = yCur < yMin[i] ? yCur : yMin[i];


         // michalak
          double yStar = 1.5;
          double yRel = y[i] / yStar;
          double yCur = y[i] < 1.0 ? y[i] + (3.0 - 2.0 * yStar) * yRel * yRel + (yStar - 2.0) * yRel * yRel * yRel : 1.0;
          yMin[i] = yCur < yMin[i] ? yCur : yMin[i];

          
          if (yMin[i] < 0)
              cout << " yMin = " << yMin[i] << " < 0 --- so strange! " << MI[i] << ' ' << mI[i] << ' ' << funval1[i] << ' ' << el_uMean[i] << ' ' << diff[i] << ' ' << mI[i] - el_uMean[i]  << ' ' << MI[i] - el_uMean[i] << "; iCell = " << iCell << "; i = " << i << endl;
      }
   } // for iPoint
}

void IndicatorBJ::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil, 
   const DenseMatrix& elfun1_mat)
{
   minValues.SetSize(mesh->GetNE());

   // loop through finite elements to compute indicator values

   mI = 1e9;
   MI = -1e9;
   yMin = 1e9;

   // get FE
   fe = fes->GetFE(iCell);
   fes->GetElementVDofs(iCell, el_vdofs);
   const int nDofs = fe->GetDof();

   // if (iCell == 1789 && myRank == 17) 
   // {
   //    cout << "elfun1_mat = ";
   //    elfun1_mat.Print(cout);
   // }

   // cout << "stencil: ";
   // stencil_num.Print(cout);

   // compute min|max values of solution for current stencil
   
   
   for (int k : stencil->cell_num)
   {
      // if (iCell == 3160) 
      // {
      //    el_uMean.Print(cout << "el_uMean for cell #" << k << ":");
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
   //    cout << "stencil = ";
   //    stencil_num.Print(cout);
   //    cout << "internal faces = ";
   //    internal_face_numbers.Print(cout);
   // }


   // run through faces again to compute values in gauss points
   for (int iFace : stencil->internal_face_numbers)
   {
      // Integration order calculation from DGTraceIntegrator
      face_el_trans = mesh->GetFaceElementTransformations(iFace);
      IntegrationPointTransformation curTrans = face_el_trans->Elem1No == iCell ? face_el_trans->Loc1 : face_el_trans->Loc2;

      int intorder = 2;//(min(face_el_trans->Elem1->OrderW(), face_el_trans->Elem2->OrderW()) +
                     //2*max(fes->GetFE(face_el_trans->Elem1No)->GetOrder(), fes->GetFE(face_el_trans->Elem2No)->GetOrder()));

      const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);
      // if (iCell == 3160) 
      // {
      //    std::cout << "tr->FaceGeom = " << face_el_trans->FaceGeom << ", intorder = " << intorder << std::endl;
      //    std::cout << "npoints = " << ir->GetNPoints() << std::endl;
      //    cout << "numFace = " << iFace << ";\n";
      // }

      updateYMin(*ir, &curTrans, elfun1_mat, iCell);

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
      //    std::cout << "tr->FaceGeom = " << face_el_trans->FaceGeom << ", intorder = " << intorder << std::endl;
      //    std::cout << "npoints = " << ir->GetNPoints() << std::endl;
      //    cout << "numFace = " << iFace << ";\n";
      // }
      updateYMin(*ir, &curTrans, elfun1_mat, iCell);
   } // for iFace

   // if (iCell == 3160) 
   // {
   //    yMin.Print(cout << "yMin = ");
   // }

   // set indicator values to the external field

   minValues[iCell] = yMin[0];//1e+9;

   for (int iEq = 0; iEq < num_equation; ++iEq)
   {
      values.GetBlock(iEq)[iCell] = yMin[iEq];
      // minValues[iCell] = yMin[iEq] < minValues[iCell] ? yMin[iEq] : minValues[iCell];
   }
};
