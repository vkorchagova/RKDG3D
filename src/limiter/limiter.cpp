#include "limiter.hpp"
  

Limiter::Limiter(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, bool _linearize, bool _haveLastHope, int _fdGroupAttribute, int _d) : 
   indicator(_ind),
   averager(_avgr),
   fes(_fes), 
   offsets(_offsets), 
   needLinearize(_linearize),
   haveLastHope(_haveLastHope),
   dim(_d) 
{
   mesh = fes->GetParMesh();
   parGridX.MakeRef(_fes, NULL);

   el_uMean.SetSize(num_equation);

   stencil = new Stencil();

   fdGroupCells.SetSize(mesh->GetNE());
   fdGroupCells = 0;

   if (_fdGroupAttribute > 0)
   {
      for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
      {
         int curAttr = mesh->GetAttribute(iCell);mesh->GetAttribute(iCell);

         if (curAttr == _fdGroupAttribute)
         {
            fdGroupCells[iCell] = 1;
         }
      }
   }
};


void Limiter::update(Vector &x)
{
   for (int sf = 0; sf < mesh->GetNSharedFaces(); sf++)
   {
      int lf = mesh->GetSharedFace(sf);
      lf2sf[lf] = sf;
   }

   xNew = x;

   parGridX.MakeRef(fes, x, 0);
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();
   averager.update(&x, &parGridX);
   averager.computeMeanValues();

   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
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

      /// Suppress slopes in defined group
      if (fdGroupCells[iCell])
      {
         averager.readElementAverageByNumber(iCell, el_uMean);
         for (int iEq = 0; iEq < num_equation; ++iEq)
         {
            for (int iDof = 0; iDof < nDofs; ++iDof)
            {
               elfun1_mat(iDof, iEq) = el_uMean(iEq);
            }
         }

         indicator.setValue(iCell, 0.0);
      }
      else // if not defined group - general limiting algorithm
      {
         // compute stencil
         getStencil(iCell);

         // check solution in cell for some troubles
         double iVal = indicator.checkDiscontinuity(iCell, stencil, elfun1_mat);

         // limit solution
         limit(iCell, iVal, nDofs, elfun1_mat);
      }
 
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
   {
      mesh->GetElementEdges(iCell,face_numbers,face_orient);
   }
   else if (dim == 3)
   {
      mesh->GetElementFaces(iCell,face_numbers,face_orient);
   }
   else
   {
      std::cout << "Unknown dimension number: " << dim << std::endl;
      return;
   }

   // current cell should be at the beginning of stencil
   stencil->cell_num.Append(iCell);

   // get table of element neighbors (OK for NC mesh)
   const Table &el_el = mesh->ElementToElementTable();
   Array<int> el2el;
   el_el.GetRow(iCell, el2el);

   for (int j = 0; j < el2el.Size(); ++j)
   {
      stencil->cell_num.Append(el2el[j]);
   }

   // if face is interior add the neighbour cell to the stencil, else remove face from stencil
   for (int iFace : face_numbers)
   {
      bool is_interior = mesh->FaceIsInterior(iFace);

      if (is_interior)
      {
         stencil->internal_face_numbers.Append(iFace);
      }
      else
      {
         int info1 = 0;
         int info2 = 0;
         int infoNC = 0;
         mesh->GetFaceInfos(iFace, &info1, &info2, &infoNC);

         if (infoNC >= 0)
         {
            stencil->internal_face_numbers.Append(iFace);
         }

         if (info2 >= 0)
         {
            std::map<int,int>::iterator mit = lf2sf.find(iFace); // find local no of face
            int sharedFaceNo = mit->second;

            stencil->shared_face_numbers.Append(sharedFaceNo);
         }
      }
   } // for iFace   
}

 
