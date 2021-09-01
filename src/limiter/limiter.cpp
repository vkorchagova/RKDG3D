#include "limiter.hpp"


Limiter::Limiter(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d) : 
   indicator(_ind),
   averager(_avgr),
   fes(_fes), 
   offsets(_offsets), 
   dim(_d) 
{
   mesh = fes->GetParMesh();
   parGridX.MakeRef(_fes, NULL);

   el_uMean.SetSize(num_equation);

   stencil = new Stencil();

   for (int sf = 0; sf < mesh->GetNSharedFaces(); sf++)
   {
      int lf = mesh->GetSharedFace(sf);
      lf2sf[lf] = sf;
   }

   cout << "create Limiter OK" << endl;
};


void Limiter::update(Vector &x)
{
   xNew = x;

   parGridX.MakeRef(x,0);
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();

   averager.update(&x, &parGridX);
   averager.computeMeanValues();


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

      indicator.checkDiscontinuity(iCell, stencil, elfun1_mat);

      for (int iEq = 0; iEq < num_equation; ++iEq)
         el_ind[iEq] = indicator.values.GetBlock(iEq)[iCell];

      limit(iCell, el_ind, elfun1_mat);

      xNew.SetSubVector(el_vdofs, el_x);

      stencil->clean();
   }

   // replace solution values to the new one
   x = xNew;
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();
}



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

 
