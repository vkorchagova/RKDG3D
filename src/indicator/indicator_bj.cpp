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
      {
         curTrans->Transform(ip, eip1);
      }
      else
      {
         eip1 = ip;
      }

      // Calculate basis functions on elements at the face
      el_shape.SetSize(fe->GetDof());
      fe->CalcShape(eip1, el_shape);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(el_shape, funval1);

      // Compute difference between average value and value in the gaussion boundary point
      averager.readElementAverageByNumber(iCell, el_uMean);
      subtract(funval1, el_uMean, diff);

      // Compute different fractions
      // for (int i = 0; i < num_equation; ++i)
      for (int i = 0; i < 1; ++i)
      {
         if (diff[i] > DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean[i])))
         {
            y[i] = (MI[i] - el_uMean[i]) / diff[i];
         }
         else if (diff[i] < - DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean[i])))
         {
            y[i] = (mI[i] - el_uMean[i]) / diff[i];
         }
         else
         {
            y[i] = 1.0;
         }

         computeCorrectionFunction(y[i],yMin[i]);
       
         if (yMin[i] < 0)
         {
            std::cout << " yMin = " << yMin[i] << " < 0 --- so strange! " << MI[i] << ' ' << mI[i] << ' ' << funval1[i] << ' ' << el_uMean[i] << ' ' << diff[i] << ' ' << mI[i] - el_uMean[i]  << ' ' << MI[i] - el_uMean[i] << "; iCell = " << iCell << "; i = " << i << std::endl;
         }
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

   // compute min|max values of solution for current stencil
   for (int k : stencil->cell_num)
   {
      averager.readElementAverageByNumber(k, el_uMean);
      for (int iSol = 0; iSol < num_equation; ++iSol)
      {
         
         mI[iSol] = el_uMean[iSol] < mI[iSol] ? el_uMean[iSol] : mI[iSol];
         MI[iSol] = el_uMean[iSol] > MI[iSol] ? el_uMean[iSol] : MI[iSol];
      }
   }

   // run through faces again to compute values in gauss points
   for (int iFace : stencil->internal_face_numbers)
   {
      // Integration order calculation from DGTraceIntegrator
      face_el_trans = mesh->GetFaceElementTransformations(iFace);
      IntegrationPointTransformation curTrans = face_el_trans->Elem1No == iCell ? face_el_trans->Loc1 : face_el_trans->Loc2;

      int intorder = 2;//(std::min(face_el_trans->Elem1->OrderW(), face_el_trans->Elem2->OrderW()) +
                      //2*std::max(fes->GetFE(face_el_trans->Elem1No)->GetOrder(), fes->GetFE(face_el_trans->Elem2No)->GetOrder()));

      const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);

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

      updateYmin(*ir, &curTrans, elfun1_mat, iCell);
   } // for iFace

   // set indicator values to the external field 
   double iVal = yMin[0];
   
   setValue(iCell, iVal);

   return iVal;
};
