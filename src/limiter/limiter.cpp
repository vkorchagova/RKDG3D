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

      // if (iCell == 1840 || iCell == 1842 || iCell == 1845)
      // {
      //   std::cout << "-- Cell #" << iCell << std::endl;
      //   std::cout << "el_x: ";
      //   el_x.Print(std::cout);
      //   std::cout << "Stencil cells: ";
      //   for (int i = 0; i < stencil->cell_num.Size(); ++i)
      //       std::cout << stencil->cell_num[i] << ' ';
      //   std::cout << std::endl;
      //   std::cout << "Stencil edges: ";
      //   for (int i = 0; i < stencil->internal_face_numbers.Size(); ++i)
      //       std::cout << stencil->internal_face_numbers[i] << ' ';
      //   std::cout << std::endl;
      // }

      // make matrix from values data
      DenseMatrix elfun1_mat(el_x.GetData(), nDofs, num_equation);

      /// Suppress slopes in defined group
      if (fdGroupCells[iCell])
      {
         averager.readElementAverageByNumber(iCell, el_uMean);
         for (int iEq = 0; iEq < num_equation; ++iEq)
            for (int iDof = 0; iDof < nDofs; ++iDof)
            {
                // if (myRank == 0 && iCell == 0 && iEq == 0) {std::cout << "   funval before lim = " << elfun1_mat(iDof, iEq) << std::endl;}
                elfun1_mat(iDof, iEq) = el_uMean(iEq);
                // if (myRank == 0 && iCell == 0 && iEq == 0) {std::cout << "   funval after lim = " << elfun1_mat(iDof, iEq) << std::endl;}
            }

         indicator.setValue(iCell, 0.0);
      }
      else // if not defined group - general limiting algorithm
      {
         // compute stencil
         getStencil(iCell);

         // check solution in cell for some troubles
         double iVal = indicator.checkDiscontinuity(iCell, stencil, elfun1_mat);

         // if ((iCell == 50 || iCell == 54) && (myRank == 0 || myRank == 3))
         // {
         //   stencil->Print(*mesh, iCell, myRank); 
         //   std::cout << "cell #" << iCell << ": iVal = "<< std::setprecision(18) << iVal << std::setprecision(6) << std::endl;
         // } 

         // if (iVal < 0.999999) iVal = 0.9; 

         limit(iCell, iVal, nDofs, elfun1_mat);

         // if (iCell == 50 || iCell == 1562)
         //   elfun1_mat.Print(std::cout << "cell #" << iCell << "; ivalue = " << iVal << "; el_fun1 = ");
      }
 
      xNew.SetSubVector(el_vdofs, el_x);

      stencil->clean();
   }

   // for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   //   if (indicator.values[iCell] < 1e-6 && myRank == 39)
   //          {
   //              std::cout << "Negative indicator value: " << indicator.values[iCell]
   //               << "; iCell = " << iCell << " in rank = " << myRank << std::endl;
   //          }

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
      std::cout << "Unknown dimension number: " << dim << std::endl;
      return;
   }

   // current cell should be at the beginning of stencil
   stencil->cell_num.Append(iCell);

   // get table of element neighbors (OK for NC mesh)
   const Table &el_el = mesh->ElementToElementTable();
   Array<int> el2el;
   el_el.GetRow(iCell, el2el);

   for (int jjj = 0; jjj < el2el.Size(); ++jjj)
      stencil->cell_num.Append(el2el[jjj]);

   // if face is interior add the neighbour cell to the stencil, else remove face from stencil
   for (int iFace : face_numbers)
   {
      // if (iCell == 1840 || iCell == 1842 || iCell == 1845)
      // {
      //   std::cout << "face #" << iFace;
      //   std::cout << " is_internal: " << mesh->FaceIsInterior(iFace);
      //   int info1 = 0;
      //   int info2 = 0;
      //   int infoNC = 0;
      //   mesh->GetFaceInfos(iFace, &info1, &info2, &infoNC);
      //   std::cout << "; info 1 = " << info1;
      //   std::cout << "; info 2 = " << info2;
      //   std::cout << "; infoNC = " << infoNC;
      //   std::cout << std::endl;
      // }
      bool is_interior = mesh->FaceIsInterior(iFace);
      //std::cout << iFace << " is interior: " << is_interior << std::endl;

      // if (myRank == 3 && iCell == 16) 
      //   {std::cout << iFace << " is interior ? " << is_interior << std::endl;}

      if (is_interior)
      {
         // face_el_trans = mesh->GetInteriorFaceTransformations(iFace);
         // int neib_cell_num = face_el_trans->Elem1No == iCell ? face_el_trans->Elem2No : face_el_trans->Elem1No;
         // if (iCell == 5369) std::cout << "my_num = " <<  (face_el_trans->Elem1No == iCell ? 1 : 2) << std::endl;

         // stencil->cell_num.Append(neib_cell_num);
         stencil->internal_face_numbers.Append(iFace);
         // stencil_size++;

         // if (myRank == 22)
         //   {
         //       std::cout << "\tinternal face_el_trans->Elem1No = " << face_el_trans->Elem1No << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << std::endl;
         //       std::cout << "\tavg neib elem = ";
         //       Vector(uMean.GetColumn(neib_cell_num), num_equation).Print(std::cout);
         //   }
      }
      else
      {
         int info1 = 0;
         int info2 = 0;
         int infoNC = 0;
         mesh->GetFaceInfos(iFace, &info1, &info2, &infoNC);

         // if (iCell == 1789 && myRank == 17)
         //   std::cout << "face no " << iFace << " with info " << info2 << std::endl;

          if (infoNC >= 0 )
          {
            stencil->internal_face_numbers.Append(iFace);
            // if (iCell == 1840 || iCell == 1842 || iCell == 1845)
            // {
            //   // for (int jjj = 0; jjj < mesh->ncmesh->GetEdgeList().masters.Size(); ++jjj)
            //   // {
            //   //   std::cout << mesh->ncmesh->GetEdgeList().masters[jjj].index << ' ';
            //   //   std::cout << mesh->ncmesh->GetEdgeList().masters[jjj].element << ' ';
            //   //   std::cout << mesh->ncmesh->GetEdgeList().masters[jjj].slaves_begin << ' ';
            //   //   std::cout << mesh->ncmesh->GetEdgeList().masters[jjj].slaves_end << ' ';
            //   //   std::cout << std::endl;
            //   // }

            //   int masterFaceID = mesh->GetNCFacesInfo()[infoNC].MasterFace;

                

            //   int slaves_begin = mesh->ncmesh->GetEdgeList().masters[masterFaceID].slaves_begin ;
            //   int slaves_end = mesh->ncmesh->GetEdgeList().masters[masterFaceID].slaves_end ;

            //   std::cout << "get num faces = " << mesh->GetNumFaces() << std::endl;
            //   // face_el_trans = mesh->GetInteriorFaceTransformations(slaves_begin);
            //   std::cout << "masterFaceID = " << masterFaceID
            //         << "; slaves beg = " << slaves_begin 
            //         << "; is_interior = " << mesh->FaceIsInterior(slaves_begin)
            //         << "; slave master = " << mesh->ncmesh->GetEdgeList().slaves[slaves_begin].master
            //         << "; slave index = " << mesh->ncmesh->GetEdgeList().slaves[slaves_begin].index
            //         << "; slave element = " << mesh->ncmesh->GetEdgeList().slaves[slaves_begin].element
            //         // << "; face_el_trans->Elem1No = " << face_el_trans->Elem1No
            //         // << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No
            //          << std::endl;
            //   // face_el_trans = mesh->GetInteriorFaceTransformations(slaves_begin);
            //   std::cout << "slaves end = " << slaves_end-1
            //         << "; is_interior = " << mesh->FaceIsInterior(slaves_end-1)
            //         << "; slave master = " << mesh->ncmesh->GetEdgeList().slaves[slaves_end-1].master
            //         << "; slave index = " << mesh->ncmesh->GetEdgeList().slaves[slaves_begin].index
            //         << "; slave element = " << mesh->ncmesh->GetEdgeList().slaves[slaves_begin].element

            //         // << "; face_el_trans->Elem1No = " << face_el_trans->Elem1No
            //         // << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No
            //          << std::endl;
            // }
         }

         if (info2 >= 0)// && iCell == 1789 && myRank == 17)
         {
            // std::cout << "found shared face no " << iFace << " with info " << info2 << std::endl;
            std::map<int,int>::iterator mit = lf2sf.find(iFace); // find local no of face
            int sharedFaceNo = mit->second;
            // if (myRank == 22)
            // std::cout << "iFace = " << iFace << "; shared face num " << sharedFaceNo << "; total num of sh faces = " << mesh->GetNSharedFaces() << std::endl;

            face_el_trans = mesh->GetSharedFaceTransformations(sharedFaceNo);
            int neibCellNo = face_el_trans->Elem2No;

            // if (myRank == 22)
            // {
            //   std::cout << "\ttr->Elem1No = " << face_el_trans->Elem1No << "; face_el_trans->Elem2No = " << face_el_trans->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << std::endl;
            //   pleaseWriteMe = true;
            // }

            // fe = fes->GetFaceNbrFE(neibCellNo - mesh->GetNE());

            // stencil->cell_num.Append(neibCellNo);
            stencil->shared_face_numbers.Append(sharedFaceNo);
         }
      }
   } // for iFace   
}

 
