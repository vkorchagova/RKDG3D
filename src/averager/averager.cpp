#include "averager.hpp"


Averager::Averager(ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d) : 
   fes(_fes), 
   offsets(_offsets), 
   dim(_d) 
{
   mesh = fes->GetParMesh();
   //parGridX->MakeRef(_fes, NULL);

   // prepare place for the average values
   fec_avg = new DG_FECollection(0, dim);
   fes_avg = new ParFiniteElementSpace(mesh, fec_avg, num_equation);
   fes_avg_component = new ParFiniteElementSpace(mesh, fec_avg);
   offsets_avg.SetSize(num_equation + 1);

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg[k] = k * fes_avg->GetNDofs();

   
   //avgs = new ParGridFunction(fes_avg, u_block_avg->GetData());
   avgs = new ParGridFunction(fes_avg);
   u_block_avg = new BlockVector(*avgs, offsets_avg);


   // prepare place for the average values
   fec_avg_extrap = new DG_FECollection(0, dim);
   fes_avg_extrap = new ParFiniteElementSpace(mesh, fec_avg_extrap, num_equation);
   offsets_avg_extrap.SetSize(num_equation + 1);

   fes_avg_extrap_component = new ParFiniteElementSpace(mesh, fec_avg_extrap);

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg_extrap[k] = k * fes_avg_extrap->GetNDofs();

   avgs_extrap = new ParGridFunction(fes_avg_extrap);
   u_block_avg_extrap = new BlockVector(*avgs_extrap, offsets_avg_extrap);
};

Averager::~Averager()
{
   delete avgs_extrap; 
   delete u_block_avg_extrap; 
   delete fes_avg_extrap_component;
   delete fes_avg_extrap; 
   delete fec_avg_extrap;
   delete avgs; 
   delete u_block_avg; 
   delete fes_avg_component;
   delete fes_avg; 
   delete fec_avg;
};


void Averager::updateSpaces()
{
   // std::cout << "before fes_avg->Update();" << std::endl;

   fes_avg->Update();
   fes_avg_component->Update();

   // std::cout << "after fes_avg->Update();" << std::endl;
   fes_avg_extrap->Update();
   fes_avg_extrap_component->Update();
}

void Averager::updateSolutions()
{
   // Update the space: recalculate the number of DOFs and construct a matrix
   // that will adjust any GridFunctions to the new mesh state.
   
   // Interpolate the solution on the new mesh by applying the transformation
   // matrix computed in the finite element space. Multiple GridFunctions could
   // be updated here.
   //parGridX->Update();
   avgs->Update();
   avgs_extrap->Update();

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg[k] = k * fes_avg->GetNDofs();

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg_extrap[k] = k * fes_avg_extrap->GetNDofs();

   u_block_avg->Update(*avgs,offsets_avg);
   u_block_avg_extrap->Update(*avgs_extrap, offsets_avg_extrap);
}

void Averager::updateFinished()
{
   // Free any transformation matrices to save memory.
   fes_avg->UpdatesFinished();
   fes_avg_component->UpdatesFinished();
   fes_avg_extrap->UpdatesFinished();
   fes_avg_extrap_component->UpdatesFinished();
}

void Averager::computeMeanValues()
{
   // std::cout << "in computeMeanValues..." << std::endl;
   
   for (int i = 0; i < num_equation; ++i)
   {
      // std::cout << "before ParGridFunction..." << std::endl;
      // std::cout << "v.Size() = " 
      //      << parGridX->Size() 
      //      << ", v_offset = " 
      //      << offsets[i]
      //      << ", f->GetVSize() = "
      //      << fes->GetVSize()
      //      << std::endl;
      ParGridFunction sol_i, avgs_i;
      // std::cout << "sol make ref before..." << std::endl;  
      sol_i.MakeRef(fes, *parGridX, offsets[i]);
      // std::cout << "sol make ref..." << std::endl; 
      avgs_i.MakeRef(fes_avg_component, *avgs, offsets_avg[i]);
      // std::cout << i << " " << sol_i.Size() << " " << avgs_i.Size() << " " << offsets[i] << " "<< offsets_avg[i] <<std::endl;
      sol_i.GetElementAverages(avgs_i);
      // std::cout << "OK" << std::endl;
   }

   avgs->ExchangeFaceNbrData();
}

void Averager::readElementAverageByNumber(const int iCell, Vector& el_uMean)
{
   if (iCell < mesh->GetNE())
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (*avgs)[iEq * mesh->GetNE() + iCell];
   else
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (avgs->FaceNbrData())[iEq  + (iCell - mesh->GetNE()) * num_equation];
}

void Averager::readElementExtrapAverageByNumber(const int iCell, Vector& el_uMean)
{
   if (iCell < mesh->GetNE())
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (*avgs_extrap)[iEq * mesh->GetNE() + iCell];
   else
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (avgs_extrap->FaceNbrData())[iEq  + (iCell - mesh->GetNE()) * num_equation];
}


void Averager::computeStencilExtrapAveragesVector(
   const Stencil* stencil  
)
{
   for (int i = 0; i < num_equation; ++i)
   {
      // std::cout << "-------------------- num eqn = " << i << std::endl;
      ParGridFunction sol_i, avgs_extrap_i;
      // std::cout << "sol make ref..." << std::endl;  
      sol_i.MakeRef(fes, *parGridX, offsets[i]);
      // std::cout << "sol avg make ref..." << std::endl; 
      avgs_extrap_i.MakeRef(fes_avg_extrap_component, *avgs_extrap, offsets_avg[i]);
      // std::cout << sol_i.Size() << " " << avgs_extrap_i.Size() << " " << offsets[i] << " "<< offsets_avg[i] <<std::endl;
      computeStencilExtrapAverages(sol_i, stencil, avgs_extrap_i);
   }

   avgs_extrap->ExchangeFaceNbrData();
}


void Averager::computeStencilExtrapAverages(
   const ParGridFunction& x,
   const Stencil* stencil,
   ParGridFunction& avgs_local
)
{
   // std::cout << "in extrap avg" << std::endl;
   MassIntegrator Mi;
   DenseMatrix loc_mass, loc_mass_shifted;
   Array<int> te_dofs, tr_dofs;
   Vector loc_avgs, loc_this;
   Vector int_psi(avgs_local.Size());

   int_psi = 0.0;

   int iCellTroubled = stencil->cell_num[0];
   // std::cout << "=== first cycle in stencil\n";

   for (int i : stencil->cell_num)
   {
      // std::cout << "i = " << i << std::endl;
      // loc_mass = shape function in gauss points * gauss weights
      Mi.AssembleElementMatrix2(
         *fes->GetFE(i), 
         *avgs_local.FESpace()->GetFE(i),
         *fes->GetElementTransformation(i), 
         loc_mass);
      // loc_mass_shifted = shape function of neighbour in gauss points of troubled element * gauss weights
      assembleShiftedElementMatrix(
         *fes->GetFE(i), 
         *fes->GetFE(iCellTroubled),
         *avgs_local.FESpace()->GetFE(i),
         *fes->GetElementTransformation(i),
         *fes->GetElementTransformation(iCellTroubled), 
         loc_mass_shifted);
      // std::cout << "after ass shift ang" << std::endl;
      fes->GetElementDofs(i, tr_dofs);
      avgs_local.FESpace()->GetElementDofs(i, te_dofs);
      x.GetSubVector(tr_dofs, loc_this);
      loc_avgs.SetSize(te_dofs.Size());
      loc_mass_shifted.Mult(loc_this, loc_avgs);
      // std::cout << "loc_avgs = ";
      // loc_avgs.Print (std::cout);
      // std::cout << "te_dofs = ";
      // te_dofs.Print (std::cout);

      avgs_local.SetSubVector(te_dofs, loc_avgs);
      // std::cout << "add elem vector: avgs = " << avgs(i) << ' '  << avgs(te_dofs[0])<< std::endl;
      loc_this = 1.0; // assume the local basis for 'this' sums to 1
      loc_mass.Mult(loc_this, loc_avgs);
      int_psi.SetSubVector(te_dofs, loc_avgs);
      // std::cout << "1st cycle: avgs / psi = " << avgs(i) << ' ' << int_psi(i) << std::endl;
   }
   // std::cout << "=== second cycle in stencil\n";
   for (int i : stencil->cell_num)
   {
      // std::cout << "before: avgs / psi = " << avgs(i) << ' ' << int_psi(i) << std::endl;
      avgs_local(i) /= int_psi(i);
      // std::cout << "after: avgs = " << avgs(i) << '\n';
   }
   // std::cout << std::endl;

   // std::cout << "======\nend extrap avg" << std::endl;
}


void Averager::assembleShiftedElementMatrix(
   const FiniteElement &trial_fe, 
   const FiniteElement &troubled_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans, 
   ElementTransformation &TransTroubled,
   DenseMatrix &elmat
)
{
   // std::cout << "in ass shift avg" << std::endl;
   int nDofsTrial = trial_fe.GetDof();
   int nDofsTroubled = troubled_fe.GetDof();
   int nDofsTest = test_fe.GetDof();
   double w;

   //#ifdef MFEM_THREAD_SAFE
   Vector trial_shape (nDofsTroubled);
   Vector te_shape (nDofsTest);
   Vector ipPhys;
   IntegrationPoint ipRef;
   // Vector* trial_shape = new Vector(nDofsTrial);
   // Vector* te_shape = new Vector(nDofsTest);
   //#endif
   elmat.SetSize(nDofsTest, nDofsTroubled);


   const int order = troubled_fe.GetOrder() + test_fe.GetOrder() + Trans.OrderW();
   const IntegrationRule *ir = &IntRules.Get(troubled_fe.GetGeomType(), order);

   elmat = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      TransTroubled.Transform(ir->IntPoint(i), ipPhys);
      Trans.TransformBack(ipPhys,ipRef);
      // std::cout << ip.x << ' ' << ip.y << std::endl;;
      trial_fe.CalcShape(ipRef, trial_shape);  // element - shifted, gauss points - troubled
      // std::cout << "trial shape = ";
      // trial_shape.Print(std::cout);
      test_fe.CalcShape(ipRef, te_shape);   // should be only 1
      // std::cout << "test shape = ";
      // te_shape.Print(std::cout);

      //Trans.SetIntPoint (&ipRef);
      w = Trans.Weight() * ip.weight;
      te_shape *= w;
      AddMultVWt(te_shape, trial_shape, elmat);
   }

   // std::cout << "elmat = ";
   //   (elmat).Print(std::cout);

   // delete te_shape;
   // delete trial_shape;
   // std::cout << "end ass shift avg" << std::endl;
}
