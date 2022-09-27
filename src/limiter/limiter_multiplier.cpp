#include "limiter_multiplier.hpp"
#include "physics.hpp"

// bool pleaseWriteMe = false;

LimiterMultiplier::LimiterMultiplier
(
   Indicator& _ind, 
   Averager& _avgr, 
   ParFiniteElementSpace*_fes, 
   const Array<int>& _offsets, 
   bool _linearize, 
   bool _haveLastHope, 
   int _fdGroupAttribute, 
   int _d
) : Limiter(_ind, _avgr, _fes, _offsets,_linearize,_haveLastHope, _fdGroupAttribute, _d) 
{
   el_shape.SetSize(num_equation);
};

void LimiterMultiplier::limit(const int iCell, const double ind_value, const double nDofs, DenseMatrix& elfun1_mat) 
{
   if (ind_value > DEFAULT_INDICATOR_CORRECTION_THRESHOLD_VALUE) 
   {
      return;
   }

   averager.readElementAverageByNumber(iCell, el_uMean);

   // replace solution to mean values
   for (int iEq = 0; iEq < num_equation; ++iEq)
   {
      for (int iDof = 0; iDof < nDofs; ++iDof)
      {
         if (needLinearize && fe->GetGeomType() != Geometry::TRIANGLE)
         {
            linearize(iCell, el_uMean, elfun1_mat);
         }

         for (int iDof = 0; iDof < nDofs; ++iDof)
         {
            elfun1_mat(iDof, iEq) = ind_value * elfun1_mat(iDof, iEq) + (1.0 - ind_value) * el_uMean[iEq];
         }
      }
   }
   
   // // Last hope limiter
   if (haveLastHope)
   {
      const IntegrationRule *irVertices = Geometries.GetVertices(fe->GetGeomType());

      for (int iPoint = 0; iPoint < irVertices->GetNPoints(); ++iPoint)
      {
         const IntegrationPoint &ip = irVertices->IntPoint(iPoint);

         el_shape.SetSize(nDofs);

         // Calculate basis functions on elements at the face
         fe->CalcShape(ip, el_shape);

         Vector funval1Vert(num_equation);

         // Interpolate elfun at the point
         elfun1_mat.MultTranspose(el_shape, funval1Vert);

         if (!StateIsPhysical(funval1Vert,dim))
         {
            averager.readElementAverageByNumber(iCell, el_uMean);

            for (int iEq = 0; iEq < num_equation; ++iEq)
            {
               for (int iDof = 0; iDof < nDofs; ++iDof) 
               { 
                  elfun1_mat(iDof, iEq) = el_uMean(iEq);
               }
            }
         }
      }
   }
};

void LimiterMultiplier::linearize(const int iCell, const Vector& el_uMean, DenseMatrix &elfun1_mat)
{
   // for gradients
   DenseMatrix dshapedr;
   DenseMatrix dF;
   Vector ip_phys(dim);

   /// Center of FE in physical space
   Vector el_center_phys;

   /// Center of FE in reference space
   IntegrationPoint el_center_ref;
   
   // get FE dofs indices
   fes->GetElementVDofs(iCell, el_vdofs);
   const int nDofs = fe->GetDof();

   // get FE
   fe = fes->GetFE(iCell);

   // set place for central gradients
   dshapedr.SetSize(nDofs, dim);
   dF.SetSize(num_equation, dim);

   // get FE transformation
   el_trans = mesh->GetElementTransformation(iCell);

   // get center of finite element
   mesh->GetElementCenter(iCell, el_center_phys);

   // move it to the element space
   el_trans->TransformBack(el_center_phys, el_center_ref);
   
   // compute gradient of el_shape functions in the center (OK for p = 1 because of 1,x,y,xy)
   fe->CalcDShape(el_center_ref, dshapedr);
   MultAtB(elfun1_mat, dshapedr, dF);

   const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), fe->GetOrder() + 1); 

   for (int iDof = 0; iDof < nDofs; ++iDof)
   {
      const IntegrationPoint ip = ir->IntPoint(iDof);
      el_trans->Transform(ip,ip_phys);

      for (int iEq = 0; iEq < num_equation; ++iEq)
      {
         elfun1_mat(iDof, iEq) = el_uMean(iEq) - \
            dF(iEq,0) * (el_center_ref.x - ip.x) - \
            dF(iEq,1) * (el_center_ref.y - ip.y) ;

         if (dim == 3)
         {
            elfun1_mat(iDof, iEq) -= dF(iEq,2) * (el_center_ref.z - ip.z) ;
         }
      }
   }
}
