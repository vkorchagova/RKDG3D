#include "indicator_bj.hpp"

// bool pleaseWriteMe = false;

IndicatorBJ::IndicatorBJ(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata) 
    : Indicator(_fes, _offsets, _d, _idata) 
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
   const ParGridFunction* uMean, 
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

      // cout << "eip in ref space = " << eip1.x << " " << eip1.y << " " << endl;
      // cout << "fe = " << fe << endl;

      // Calculate basis functions on elements at the face
      el_shape.SetSize(fe->GetDof());
      fe->CalcShape(eip1, el_shape);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(el_shape, funval1);
      // cout << "funval1: ";
      // funval1.Print(cout);
      // cout << "uMean[iCell]: ";
      // uMean[iCell].Print(cout);

      // if (iCell == 1789 && myRank == 17) 
      // {
      //    // cout << "eip in ref space = " << eip1.x << " " << eip1.y << " " << endl;
      //    // cout << "shape = ";
      //    // el_shape.Print(cout);
      //    // cout << "elfun1_mat = ";
      //    // elfun1_mat.Print(cout);
      //    cout << "funval1_face = ";
      //    funval1.Print(cout);
      // }


      //uMean.GetCol(iCell,uMean_el);
      // for (int iii = 0; iii < funval1.Size(); ++iii)
      //    if (funval1[iii] != funval1[iii])
      //    {
      //       cout << "Find NaN for funval1[" << iii << "], processor " << myRank;
      //       cout << "; iCell = " << iCell << "; face_el_trans elem 1 no = "<< face_el_trans->Elem1No << "; face_el_trans elem 2 no =" << face_el_trans->Elem2No << endl;
      //    }
      // if (iCell == 1789 && myRank == 17) 
      // {
      //    cout << "uMean = ";
      //    Vector(uMean.GetColumn(iCell), num_equation).Print(cout);
      // }

      readElementAverageByNumber(iCell, mesh, uMean, el_uMean);
      subtract(funval1, el_uMean, diff);
      // if (iCell == 5369)
      // {
      //    cout << "\tgp in ref space = " << eip1.x << " " << eip1.y << " " ;
      //    cout << "funVal_faceGP = ";
      //    funval1.Print(cout);
      // }

      // Compute different fractions
      for (int i = 0; i < num_equation; ++i)
      {
         if (diff[i] > 1e-6)
              y[i] = (MI[i] - el_uMean[i]) / diff[i];
          else if (diff[i] < -1e-6)
              y[i] = (mI[i] - el_uMean[i]) / diff[i];
          else
              y[i] = 1.0;
         
         // // classical BJ
         yMin[i] = y[i] < yMin[i] ? y[i] : yMin[i];
         yMin[i] = yMin[i] > 1.0 ? 1.0 : yMin[i];

         // // venkatakrishnan
         // double yCur = (y[i]*y[i] + 2.0 * y[i]) / (y[i] * y[i] + y[i] + 2.0);
         // yMin[i] = yCur < yMin[i] ? yCur : yMin[i];


         // michalak
         // double yStar = 1.5;
         // double yRel = y[i] / yStar;
         // double yCur = y[i] < 1.0 ? y[i] + (3.0 - 2.0 * yStar) * yRel * yRel + (yStar - 2.0) * yRel * yRel * yRel : 1.0;
         // yMin[i] = yCur < yMin[i] ? yCur : yMin[i];

          
          if (yMin[i] < 0)
              cout << " yMin = " << yMin[i] << " < 0 --- so strange! " << MI[i] << ' ' << mI[i] << ' ' << funval1[i] << ' ' << el_uMean[i] << ' ' << diff[i] << ' ' << mI[i] - el_uMean[i]  << ' ' << MI[i] - el_uMean[i] << "; iCell = " << iCell << "; i = " << i << endl;
      }
   } // for iPoint
}

void IndicatorBJ::checkDiscontinuity(
   const int iCell, 
   const Stencil* stencil, 
   const ParGridFunction* uMean, 
   const DenseMatrix& elfun1_mat)
{
   Vector uMean_el(num_equation);

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
         for (int iSol = 0; iSol < num_equation; ++iSol)
         {
            readElementAverageByNumber(k, mesh, uMean, el_uMean);
            mI[iSol] = el_uMean[iSol] < mI[iSol] ? el_uMean[iSol] : mI[iSol];
            MI[iSol] = el_uMean[iSol] > MI[iSol] ? el_uMean[iSol] : MI[iSol];
         }

      // if (iCell == 1789 && myRank == 17) 
      // {
      //    cout << "mI = ";
      //    mI.Print(cout);
      //    cout << "MI = ";
      //    MI.Print(cout);
      // }


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

         int intorder = (min(face_el_trans->Elem1->OrderW(), face_el_trans->Elem2->OrderW()) +
                        2*max(fes->GetFE(face_el_trans->Elem1No)->GetOrder(), fes->GetFE(face_el_trans->Elem2No)->GetOrder()));

         const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);
         // if (iCell == 1789 && myRank == 17) 
         // {
         //    std::cout << "tr->FaceGeom = " << face_el_trans->FaceGeom << ", intorder = " << intorder << std::endl;
         //    std::cout << "npoints = " << ir->GetNPoints() << std::endl;
         //    cout << "numFace = " << iFace << ";\n";
         // }

         updateYMin(*ir, &curTrans, elfun1_mat, uMean, iCell);

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
         updateYMin(*ir, &curTrans, elfun1_mat, uMean, iCell);
      } // for iFace

      // set indicator values

      for (int iEq = 0; iEq < num_equation; ++iEq)
         values.GetBlock(iEq)[iCell] = yMin[iEq];
};
