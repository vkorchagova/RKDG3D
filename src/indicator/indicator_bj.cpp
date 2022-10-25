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
   // Vector diff(num_equation);
   // Vector y(num_equation);

   double diff;
   double y;

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
      averager.readElementAverageComponentByNumber(iCell, iSolRoot, el_uMean);
      diff = funval1[iSolRoot] - el_uMean;
      
      if (diff > DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean)))
      {
         y = (MI - el_uMean ) / diff;
      }
      else if (diff< - DEFAULT_BJ_DIFF_MAX_PERCENT * std::max(1.0, fabs(el_uMean )))
      {
         y = (mI - el_uMean ) / diff;
      }
      else
      {
         y = 1.0;
      }

      computeCorrectionFunction(y,yMin);
      
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
      averager.readElementAverageComponentByNumber(k, iSolRoot, el_uMean);
      mI = el_uMean < mI ? el_uMean : mI;
      MI = el_uMean > MI ? el_uMean : MI;
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

      int intorder = face_el_trans->Elem1->OrderW() + 2*fes->GetFE(face_el_trans->Elem1No)->GetOrder();

      const IntegrationRule *ir = &IntRules.Get(face_el_trans->FaceGeom, intorder);

      updateYmin(*ir, &curTrans, elfun1_mat, iCell);
   } // for iFace

   // set indicator values to the external field 
   double iVal = yMin;

   setValue(iCell, iVal);

   return iVal;
};
