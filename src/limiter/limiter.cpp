#include "limiter.hpp"

void Limiter::GetStencilAverages(
   const GridFunction& x,
   const Array<int>& stencil_num,
   GridFunction& avgs
) const
{
   MassIntegrator Mi;
   DenseMatrix loc_mass;
   Array<int> te_dofs, tr_dofs;
   Vector loc_avgs, loc_this;
   Vector int_psi(avgs.Size());

   avgs = 0.0;
   int_psi = 0.0;

   int iCellTroubled = stencil_num[0];

   for (int i : stencil_num)
   {
      // loc_mass = shape function in gauss points * gauss weights

      AssembleShiftedElementMatrix(
         *fes->GetFE(i), 
         *fes->GetFE(iCellTroubled),
         *avgs.FESpace()->GetFE(i),
         *fes->GetElementTransformation(i), 
         loc_mass);

      fes->GetElementDofs(i, tr_dofs);
      avgs.FESpace()->GetElementDofs(i, te_dofs);
      x.GetSubVector(tr_dofs, loc_this);
      loc_avgs.SetSize(te_dofs.Size());
      loc_mass.Mult(loc_this, loc_avgs);

      avgs.AddElementVector(te_dofs, loc_avgs);
      loc_this = 1.0; // assume the local basis for 'this' sums to 1
      loc_mass.Mult(loc_this, loc_avgs);
      int_psi.AddElementVector(te_dofs, loc_avgs);
   }
   for (int i = 0; i < avgs.Size(); i++)
   {
      avgs(i) /= int_psi(i);
   }
}


void Limiter::AssembleShiftedElementMatrix(
   const FiniteElement &trial_fe, 
   const FiniteElement &troubled_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans, 
   DenseMatrix &elmat) const
{
   int nDofsTrial = trial_fe.GetDof();
   int nDofsTroubled = troubled_fe.GetDof();
   int nDofsTest = test_fe.GetDof();
   double w;

//#ifdef MFEM_THREAD_SAFE
   Vector trial_shape, te_shape;
//#endif
   elmat.SetSize(nDofsTest, nDofsTrial);
   trial_shape.SetSize(nDofsTrial);
   te_shape.SetSize(nDofsTroubled);

   const int order = trial_fe.GetOrder() + test_fe.GetOrder() + Trans.OrderW();
   const IntegrationRule *ir = &IntRules.Get(trial_fe.GetGeomType(), order);

   elmat = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      troubled_fe.CalcShape(ip, trial_shape);  // element - troubled, gauss points - shifted
      test_fe.CalcShape(ip, te_shape);   // should be only 1

      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      te_shape *= w;
      AddMultVWt(te_shape, trial_shape, elmat);
   }
}

 
