#include "limiter.hpp"


Limiter::Limiter(Indicator& _ind, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d) : 
   indicator(_ind),
   fes(_fes), 
   offsets(_offsets), 
   dim(_d) 
{
   mesh = fes->GetParMesh();
   parGridX.MakeRef(_fes, NULL);

   // prepare place for the average values
   fec_avg = new DG_FECollection(0, dim);
   fes_avg = new ParFiniteElementSpace(mesh, fec_avg, num_equation);
   offsets_avg.SetSize(num_equation + 1);

   for (int k = 0; k <= num_equation; k++) 
      offsets_avg[k] = k * fes_avg->GetNDofs();

   u_block_avg = new BlockVector(offsets_avg);
   avgs = new ParGridFunction(fes_avg, u_block_avg->GetData());

   el_uMean.SetSize(num_equation);

   stencil = new Stencil();

   for (int sf = 0; sf < mesh->GetNSharedFaces(); sf++)
   {
      int lf = mesh->GetSharedFace(sf);
      lf2sf[lf] = sf;
   }
};


void Limiter::update(Vector &x)
{
   xNew = x;

   parGridX.MakeRef(x,0);
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();

   computeMeanValues();

   Vector el_ind(num_equation);

   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
      // compute stencil
      getStencil(iCell);

      // get FE
      fe = fes->GetFE(iCell);

      // get FE dofs indices
      fes->GetElementVDofs(iCell, el_vdofs);
      const int nDofs = fe->GetDof();

      //get FE transformation
      el_trans = mesh->GetElementTransformation(iCell);
      //el_x.SetSize(nDofs);

      // read sol values for defined el_vdofs
      x.GetSubVector(el_vdofs, el_x);

      // make matrix from values data
      DenseMatrix elfun1_mat(el_x.GetData(), nDofs, num_equation);

      indicator.checkDiscontinuity(iCell, stencil, avgs, elfun1_mat);

      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_ind[iEq] = indicator.values.GetBlock(iEq)[iCell];

      limit(iCell, el_ind, elfun1_mat);

      xNew.SetSubVector(el_vdofs, el_x);

      cleanStencil();
   }

   // replace solution values to the new one
   x = xNew;
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();
}

void Limiter::computeMeanValues()
{
   ParFiniteElementSpace fes_avg_component(mesh, fec_avg);
   for (int i = 0; i < num_equation; ++i)
   {
      ParGridFunction sol_i, avgs_i;
      // cout << "sol make ref..." << endl;  
      sol_i.MakeRef(fes, parGridX, offsets[i]);
      // cout << "sol make ref..." << endl; 
      avgs_i.MakeRef(&fes_avg_component, *avgs, offsets_avg[i]);
      // cout << i << " " << sol_i.Size() << " " << avgs_i.Size() << " " << offsets[i] << " "<< offsets_avg[i] <<endl;
      sol_i.GetElementAverages(avgs_i);
      // cout << avgs_i[6661] << endl;
   }

   avgs->ExchangeFaceNbrData();
}

void readElementAverageByNumber(const int iCell, const ParMesh* mesh, const ParGridFunction* avgs, Vector& el_uMean)
{
   if (iCell < mesh->GetNE())
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (*avgs)[iEq * mesh->GetNE() + iCell];
   else
      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_uMean[iEq] = (avgs->FaceNbrData())[iEq  + (iCell - mesh->GetNE()) * num_equation];
}


// void getStencilAverages(
//    const GridFunction& x,
//    const Array<int>& cell_num,
//    GridFunction& avgs
// )
// {
//    MassIntegrator Mi;
//    DenseMatrix loc_mass;
//    Array<int> te_dofs, tr_dofs;
//    Vector loc_avgs, loc_this;
//    Vector int_psi(avgs.Size());

//    avgs = 0.0;
//    int_psi = 0.0;

//    int iCellTroubled = cell_num[0];

//    for (int i : cell_num)
//    {
//       // loc_mass = shape function in gauss points * gauss weights

//       AssembleShiftedElementMatrix(
//          *fes->GetFE(i), 
//          *fes->GetFE(iCellTroubled),
//          *avgs.FESpace()->GetFE(i),
//          *fes->GetElementTransformation(i), 
//          loc_mass);

//       fes->GetElementDofs(i, tr_dofs);
//       avgs.FESpace()->GetElementDofs(i, te_dofs);
//       x.GetSubVector(tr_dofs, loc_this);
//       loc_avgs.SetSize(te_dofs.Size());
//       loc_mass.Mult(loc_this, loc_avgs);

//       avgs.AddElementVector(te_dofs, loc_avgs);
//       loc_this = 1.0; // assume the local basis for 'this' sums to 1
//       loc_mass.Mult(loc_this, loc_avgs);
//       int_psi.AddElementVector(te_dofs, loc_avgs);
//    }
//    for (int i = 0; i < avgs.Size(); i++)
//    {
//       avgs(i) /= int_psi(i);
//    }
// }


// void assembleShiftedElementMatrix(
//    const FiniteElement &trial_fe, 
//    const FiniteElement &troubled_fe,
//    const FiniteElement &test_fe,
//    ElementTransformation &Trans, 
//    DenseMatrix &elmat) const
// {
//    int nDofsTrial = trial_fe.GetDof();
//    int nDofsTroubled = troubled_fe.GetDof();
//    int nDofsTest = test_fe.GetDof();
//    double w;

//    //#ifdef MFEM_THREAD_SAFE
//    Vector trial_shape, te_shape;
//    //#endif
//    elmat.SetSize(nDofsTest, nDofsTrial);
//    trial_shape.SetSize(nDofsTrial);
//    te_shape.SetSize(nDofsTroubled);

//    const int order = trial_fe.GetOrder() + test_fe.GetOrder() + Trans.OrderW();
//    const IntegrationRule *ir = &IntRules.Get(trial_fe.GetGeomType(), order);

//    elmat = 0.0;
//    for (int i = 0; i < ir->GetNPoints(); i++)
//    {
//       const IntegrationPoint &ip = ir->IntPoint(i);
//       troubled_fe.CalcShape(ip, trial_shape);  // element - troubled, gauss points - shifted
//       test_fe.CalcShape(ip, te_shape);   // should be only 1

//       Trans.SetIntPoint (&ip);
//       w = Trans.Weight() * ip.weight;
//       te_shape *= w;
//       AddMultVWt(te_shape, trial_shape, elmat);
//    }
// }

void Limiter::getStencil(const int iCell)
{
   Array<int> face_numbers;
   Array<int> face_orient;

   // find faces (3D) or edges (2D) for element
   if (dim == 2)
      mesh->GetElementEdges(iCell,face_numbers,face_orient);
   else if (dim == 3)
      mesh->GetElementFaces(iCell,face_numbers,face_orient);
   else
   {
      cout << "Unknown dimension number: " << dim << endl;
      return;
   }

   // current cell should be at the beginning of stencil
   stencil->cell_num.Append(iCell);

   // if face is interior add the neighbour cell to the stencil, else remove face from stencil
   for (int iFace : face_numbers)
   {
      bool is_interior = mesh->FaceIsInterior(iFace);
      //std::cout << iFace << " is interior: " << is_interior << std::endl;

      // if (myRank == 3 && iCell == 16) 
      //    {std::cout << iFace << " is interior ? " << is_interior << std::endl;}

      if (is_interior)
      {
         face_el_trans = mesh->GetInteriorFaceTransformations(iFace);
         int neib_cell_num = face_el_trans->Elem1No == iCell ? face_el_trans->Elem2No : face_el_trans->Elem1No;
         // if (iCell == 5369) cout << "my_num = " <<  (face_el_trans->Elem1No == iCell ? 1 : 2) << endl;

         stencil->cell_num.Append(neib_cell_num);
         stencil->internal_face_numbers.Append(iFace);
         // stencil_size++;

         // if (myRank == 22)
         //    {
         //       cout << "\tinternal face_el_trans->Elem1No = " << face_el_trans->Elem1No << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << endl;
         //       cout << "\tavg neib elem = ";
         //       Vector(uMean.GetColumn(neib_cell_num), num_equation).Print(cout);
         //    }
      }
      else
      {
         int info1 = 0;
         int info2 = 0;
         mesh->GetFaceInfos(iFace, &info1, &info2);

         // if (iCell == 1789 && myRank == 17)
         //    cout << "face no " << iFace << " with info " << info2 << endl;

         if (info2 >= 0)// && iCell == 1789 && myRank == 17)
         {
            // cout << "found shared face no " << iFace << " with info " << info2 << endl;
            std::map<int,int>::iterator mit = lf2sf.find(iFace); // find local no of face
            int sharedFaceNo = mit->second;
            // if (myRank == 22)
            // cout << "iFace = " << iFace << "; shared face num " << sharedFaceNo << "; total num of sh faces = " << mesh->GetNSharedFaces() << endl;

            face_el_trans = mesh->GetSharedFaceTransformations(sharedFaceNo);
            int neibCellNo = face_el_trans->Elem2No;

            // if (myRank == 22)
            // {
            //    cout << "\ttr->Elem1No = " << face_el_trans->Elem1No << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << endl;
            //    pleaseWriteMe = true;
            // }

            // fe = fes->GetFaceNbrFE(neibCellNo - mesh->GetNE());



            stencil->cell_num.Append(neibCellNo);
            stencil->shared_face_numbers.Append(sharedFaceNo);
         }
      }
   } // for iFace
}

void Limiter::cleanStencil()
{
   stencil->internal_face_numbers.DeleteAll();
   stencil->shared_face_numbers.DeleteAll();
   stencil->cell_num.DeleteAll();
}
 
