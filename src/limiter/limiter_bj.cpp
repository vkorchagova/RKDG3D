#include "limiter_bj.hpp"
#include "physics.hpp"

bool pleaseWriteMe = false;

LimiterBJ::LimiterBJ(ParFiniteElementSpace*_fes, const Array<int>& _offsets, int _d) : Limiter(_fes, _offsets, _d) 
{
   mesh = fes->GetParMesh();
   uMean.SetSize(num_equation, mesh->GetNE() + mesh->GetNSharedFaces());
   stencil.SetSize(mesh->GetNE());

   mI.SetSize(num_equation);
   MI.SetSize(num_equation);

   // for (int i = 0; i < mesh->GetNE(); ++i)
   //    uMean[i].SetSize(num_equation);

   // prepare place for the average values
   fec_avg = new DG_FECollection(0, dim);
   fes_avg = new ParFiniteElementSpace(mesh, fec_avg, num_equation);
   offsets_avg.SetSize(num_equation + 1);
   
   for (int k = 0; k <= num_equation; k++) 
      offsets_avg[k] = k * fes_avg->GetNDofs();

   u_block_avg = new BlockVector(offsets_avg);
   avgs = new ParGridFunction(fes_avg, u_block_avg->GetData());

   // compute stencils
   //Init(mesh);
};

void LimiterBJ::Init(const ParMesh* mesh)
{
   // for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   // {
   //    // find faces (3D) or edges (2D) for element
   //    if (dim == 2)
   //       mesh->GetElementEdges(iCell,face_numbers,face_orient);
   //    else if (dim == 3)
   //       mesh->GetElementFaces(iCell,face_numbers,face_orient);
   //    else
   //    {
   //       cout << "Unknown dimension number: " << dim << endl;
   //       return;
   //    }

   //    // current cell should be at the beginning of stencil
   //    stencil[iCell].Append(iCell);

   //    // run throught faces and find cells for stencil
   //    for (int iFace : face_numbers)
   //    {
   //       bool is_interior = mesh->FaceIsInterior(iFace);
   //       //std::cout << iFace << " is interior: " << is_interior << std::endl;

   //       // if (myRank == 3 && iCell == 16) 
   //       //    {std::cout << iFace << " is interior ? " << is_interior << std::endl;}

   //       // if face is interior add the neighbor cell to the stencil
   //       if (is_interior)
   //       {
   //          tr = mesh->GetInteriorFaceTransformations(iFace);
   //          int neib_cell_num = tr->Elem1No == iCell ? tr->Elem2No : tr->Elem1No;
   //          stencil[iCell].Append(neib_cell_num);

            
   //          internal_face_numbers.Append(iFace);
   //          // stencil_size++;
   //       }
   //       else // if face is shared find the number of shared face
   //       {
   //          int info1 = 0;
   //          int info2 = 0;
   //          mesh->GetFaceInfos(iFace, &info1, &info2);

   //          if (info2 >= 0 && iCell == 1789 && myRank == 17)
   //          {
   //             cout << "found shared face no " << iFace << " with info " << info2 << endl;
   //             std::map<int,int>::iterator mit = lf2sf.find(iFace); // find local no of face
   //             int sharedFaceNo = mit->second;
   //             cout << "\tshared face num " << sharedFaceNo << "; total num of sh faces = " << mesh->GetNSharedFaces() << endl;
   //             tr = mesh->GetSharedFaceTransformations(sharedFaceNo);
   //             int neibCellNo = -1 - tr->Elem2No;
   //             cout << "tr->Elem1No = " << tr->Elem1No << endl;
   //             cout << "tr->Elem2No = " << tr->Elem2No << endl;
   //             cout << "rho_avg = " << uMean(0, tr->Elem2No) << endl;
   //             // cout << "mesh ne = " << mesh->GetNE() << endl;
   //          }

   //          //tr = mesh->GetSharedFaceTransformations(iFace);
   //       }
   //    } // for iFace

   // } // for iCell

}

void LimiterBJ::linearize(const int iCell, const DenseMatrix & uMean, DenseMatrix &elfun1_mat)
{
   // for gradients
   DenseMatrix dshapedr;
   DenseMatrix dF;
   Vector ip_phys(dim);
   
   // get FE dofs indices
   fes->GetElementVDofs(iCell, vdofs);
   const int nDofs = fe->GetDof();
   // get FE
   fe = fes->GetFE(iCell);

   // set place for central gradients
   dshapedr.SetSize(nDofs, dim);
   dF.SetSize(num_equation, dim);


   // get FE transformation
   el_trans = mesh->GetElementTransformation(iCell);

   // std::cout << "elfun1_mat = ";
   // elfun1_mat.Print(std::cout);

   // get center of finite element
   mesh->GetElementCenter(iCell, el_center_phys);
   // if (iCell == 6950) el_center_phys.Print(std::cout);

   // move it to the element space
   el_trans->TransformBack(el_center_phys, el_center_ref);
   // if (iCell == 6950) std::cout << "el_center_ref = " << el_center_ref.index << " " << el_center_ref.x << " " << el_center_ref.y << " " << el_center_ref.z << std::endl;

   // compute gradient of shape functions in the center (OK for p = 1 because of 1,x,y,xy)
   fe->CalcDShape(el_center_ref, dshapedr);
   MultAtB(elfun1_mat, dshapedr, dF);

   const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), fe->GetOrder() + 1);

   // if (iCell == 6950)
   // {
   //    cout << "etype = " << fe->GetGeomType() << "; order = " << fe->GetOrder() << endl;
   // }
   

   for (int iDof = 0; iDof < nDofs; ++iDof)
   {
      const IntegrationPoint ip = ir->IntPoint(iDof);
      el_trans->Transform(ip,ip_phys);

      // if (iCell == 6950)
      // {
      //    cout << "point = " << ip.x << ' ' << ip.y << "; phys = " << ip_phys[0] << ' ' << ip_phys[1] << endl;
      // }
   
      for (int iEq = 0; iEq < num_equation; ++iEq)
      {
         // if (iCell == 6950) cout << "\t old val = " << elfun1_mat(iDof, iEq) << endl;
         // if (iCell == 6950) cout << "\t dF = " << dF(iEq,0) << ' ' << dF(iEq,1) << endl;

         elfun1_mat(iDof, iEq) = uMean(iEq,iCell) - \
            dF(iEq,0) * (el_center_ref.x - ip.x) - \
            dF(iEq,1) * (el_center_ref.y - ip.y) ;

         if (dim == 3)
            elfun1_mat(iDof, iEq) -= dF(iEq,2) * (el_center_ref.z - ip.z) ;

            // dF(iEq,0) * (el_center_phys[0] - ip_phys[0]) - \
            dF(iEq,1) * (el_center_phys[1] - ip_phys[1]) ;
         // if (iCell == 6950) cout << "\t new val = " << elfun1_mat(iDof, iEq) << endl;
      }
   }
}



void LimiterBJ::updateYMin(const IntegrationRule& ir, IntegrationPointTransformation* curTrans, const DenseMatrix& elfun1_mat, const int iCell, Vector &yMin) 
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
      shape.SetSize(fe->GetDof());
      fe->CalcShape(eip1, shape);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(shape, funval1);
      // cout << "funval1: ";
      // funval1.Print(cout);
      // cout << "uMean[iCell]: ";
      // uMean[iCell].Print(cout);

      // if (iCell == 1789 && myRank == 17) 
      // {
      //    // cout << "eip in ref space = " << eip1.x << " " << eip1.y << " " << endl;
      //    // cout << "shape = ";
      //    // shape.Print(cout);
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
      //       cout << "; iCell = " << iCell << "; tr elem 1 no = "<< tr->Elem1No << "; tr elem 2 no =" << tr->Elem2No << endl;
      //    }
      // if (iCell == 1789 && myRank == 17) 
      // {
      //    cout << "uMean = ";
      //    Vector(uMean.GetColumn(iCell), num_equation).Print(cout);
      // }

      if (primitive)
         TransformConservativeToPrimitive(funval1);

      subtract(funval1, Vector(uMean.GetColumn(iCell), num_equation), diff);
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
              y[i] = (MI[i] - uMean(i,iCell)) / diff[i];
          else if (diff[i] < -1e-6)
              y[i] = (mI[i] - uMean(i,iCell)) / diff[i];
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
              cout << "yMin < 0 --- so strange!" << endl;
      }
   } // for iPoint
}

void LimiterBJ::limit(Vector &x) 
{
   mesh = fes->GetParMesh();
   parGridX.MakeRef(x,0);
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData(); 

   xNew = x;

   Array<int> face_numbers;
   Array<int> face_orient;

   stencil_num.Reserve(stencil_max_size);

   Vector yMin(num_equation);
   Vector uMean_el(num_equation);

   Array<int> internal_face_numbers;
   Array<int> shared_face_numbers;


   int stencil_size = 0;

   double ref = 0;

   // cout << "======\n";

   // get averages via grid fun
   ParFiniteElementSpace fes_avg_rho(mesh, fec_avg);
   for (int i = 0; i < num_equation; ++i)
   {
      ParGridFunction sol_i, avgs_i;
      // cout << "sol make ref..." << endl;  
      sol_i.MakeRef(fes, parGridX, offsets[i]);
      // cout << "sol make ref..." << endl; 
      avgs_i.MakeRef(&fes_avg_rho, *avgs, offsets_avg[i]);
      // cout << i << " " << sol_i.Size() << " " << avgs_i.Size() << " " << offsets[i] << " "<< offsets_avg[i] <<endl;
      sol_i.GetElementAverages(avgs_i);
      // cout << avgs_i[6661] << endl;
   }

   avgs->ExchangeFaceNbrData();

   // loop through finite elements to compute mean values
   // here will be only troubled cells + their neighbors after indicator
   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
      
      fe = fes->GetFE(iCell);
      //computeElementAverage(x, fe, iCell, uMean_el);
      for (int iEq = 0; iEq < num_equation; ++iEq)
         uMean_el[iEq] = (*avgs)[iEq * mesh->GetNE() + iCell];

      uMean.SetCol(iCell, uMean_el);

      // if (iCell == 1789 && myRank == 17) 
      // {
      //    cout << "uMean = ";
      //    uMean_el.Print(cout);
      // }
   }

   // loop through finite elements to compute pure gradients in centers
   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
      // get FE
      fe = fes->GetFE(iCell);

      // get FE dofs indices
      fes->GetElementVDofs(iCell, vdofs);
      const int nDofs = fe->GetDof();

      //get FE transformation
      el_trans = mesh->GetElementTransformation(iCell);
      //el_x.SetSize(nDofs);

      // read sol values for defined vdofs
      x.GetSubVector(vdofs, el_x);

      // make matrix from values data
      DenseMatrix elfun1_mat(el_x.GetData(), nDofs, num_equation);
      // std::cout << "elfun1_mat = ";
      // elfun1_mat.Print(std::cout);

      // get center of finite element
      mesh->GetElementCenter(iCell, el_center_phys);
      //if (iCell == 0) {std::cout << "el_center_phys = "; el_center_phys.Print(std::cout);}

      // move it to the element space
      el_trans->TransformBack(el_center_phys, el_center_ref);
      
      // if (iCell == 6661) 
      // {
      //    cout << "uMean_el" << endl;
      //    uMean_el.Print(cout);

      //    // cout << "dShape" << endl;
      //    // dshapedr.Print(cout);

      //    // cout << "elfun1_mat" << endl;
      //    // elfun1_mat.Print(cout);

      //    // cout << "dF" << endl;
      //    // dF.Print(cout);
      // }

      // // find gauss points and loop through them
      // const IntegrationRule *ir = &IntRules.Get(el_trans->ElementType, fe->GetOrder());

      // for (int iEq = 0; iEq < num_equation; ++iEq)
      // {
      //    for (int iDof = 0; iDof < nDofs; ++iDof)
      //    {
      //       const IntegrationPoint ip = ir->IntPoint(iDof);
      //       el_trans->Transform(ip,ip_phys);
            
      //       // if (iCell == 6661) 
      //       // {
      //       //    // cout << "===\nip = " << ip_phys[0] << " " << ip_phys[1] << endl;
      //       //    // cout << "iDof iEqn = " << iDof << " " << iEq << endl;
      //       //    cout << "elfun1_mat(" << iDof << "," << iEq <<") = " << elfun1_mat(iDof, iEq) << endl;
      //       // }

      //       elfun1_mat(iDof, iEq) = uMean(iEq,iCell) + \
      //          dF(iEq,0) * (el_center_phys[0] - ip_phys[0]) + \
      //          dF(iEq,1) * (el_center_phys[1] - ip_phys[1]);

      //       // if (iCell == 6661) 
      //       // {
      //       //    cout << "elfun1_mat(" << iDof << "," << iEq <<") new = " << elfun1_mat(iDof, iEq) << endl;

      //       // }
      //    }
         
      // }

      // x.SetSubVector(vdofs, el_x);



   } // for iCell

   // cout << "mean values are computed" << endl;
   // for (auto el_mean : uMean)
   //    el_mean.Print(cout);

   // create mapping "local-to-shared-face"
   std::map<int,int> lf2sf;
   for (int sf = 0; sf < mesh->GetNSharedFaces(); sf++)
   {
      int lf = mesh->GetSharedFace(sf);
      lf2sf[lf] = sf;
   }

   // if (myRank == 22) uMean.Print(cout);

   // loop through finite elements to compute limit values
   // here will be only troubled cells after indicator
   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
      mI = 1e9;
      MI = -1e9;
      yMin = 1e9;

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
      stencil_num.Append(iCell);

      // if (myRank == 22) cout << "ICELL = " << iCell << endl;

      // if face is interior add the neighbour cell to the stencil, else remove face from stencil
      for (int iFace : face_numbers)
      {
         bool is_interior = mesh->FaceIsInterior(iFace);
         //std::cout << iFace << " is interior: " << is_interior << std::endl;

         // if (myRank == 3 && iCell == 16) 
         //    {std::cout << iFace << " is interior ? " << is_interior << std::endl;}

         if (is_interior)
         {
            tr = mesh->GetInteriorFaceTransformations(iFace);
            int neib_cell_num = tr->Elem1No == iCell ? tr->Elem2No : tr->Elem1No;
            // if (iCell == 5369) cout << "my_num = " <<  (tr->Elem1No == iCell ? 1 : 2) << endl;

            stencil_num.Append(neib_cell_num);
            internal_face_numbers.Append(iFace);
            // stencil_size++;

            // if (myRank == 22)
            //    {
            //       cout << "\tinternal tr->Elem1No = " << tr->Elem1No << "; tr->Elem2No = " << tr->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << endl;
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

               tr = mesh->GetSharedFaceTransformations(sharedFaceNo);
               int neibCellNo = tr->Elem2No;

               // if (myRank == 22)
               // {
               //    cout << "\ttr->Elem1No = " << tr->Elem1No << "; tr->Elem2No = " << tr->Elem2No << "; mesh->GetNE() = " << mesh->GetNE() << endl;
               //    pleaseWriteMe = true;
               // }

               // fe = fes->GetFaceNbrFE(neibCellNo - mesh->GetNE());

               // computeElementAverage(x, fe, sharedFaceNo, uMean_el);
               int nNeibElems = avgs->FaceNbrData().Size();

               for (int iEq = 0; iEq < num_equation; ++iEq)
                  uMean_el[iEq] = (avgs->FaceNbrData())[iEq  + (neibCellNo - mesh->GetNE()) * num_equation];

               pleaseWriteMe = false;

               // if (myRank == 22)
               // {
               //    cout << "\tavg neib elem = ";
               //    uMean_el.Print(cout);
               // }
               uMean.SetCol(neibCellNo, uMean_el);

               stencil_num.Append(neibCellNo);
               shared_face_numbers.Append(sharedFaceNo);
            }
         }
      } // for iFace

      // get FE dofs indices
      //fes->GetElementVDofs(iCell, vdofs);

      // get FE again
      fe = fes->GetFE(iCell);
      fes->GetElementVDofs(iCell, vdofs);
      const int nDofs = fe->GetDof();

      // read sol values for defined vdofs
      x.GetSubVector(vdofs, el_x);

      // make matrix from values data
      DenseMatrix elfun1_mat(el_x.GetData(), nDofs, num_equation);

      // if (iCell == 1789 && myRank == 17) 
      // {
      //    cout << "elfun1_mat = ";
      //    elfun1_mat.Print(cout);
      // }

      // cout << "stencil: ";
      // stencil_num.Print(cout);

      // compute min|max values of solution for current stencil
      
      for (int iSol = 0; iSol < num_equation; ++iSol)
         for (int k : stencil_num)
         {
            mI[iSol] = uMean(iSol,k) < mI[iSol] ? uMean(iSol,k) : mI[iSol];
            MI[iSol] = uMean(iSol,k) > MI[iSol] ? uMean(iSol,k) : MI[iSol];
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
      for (int iFace : internal_face_numbers)
      {
         // Integration order calculation from DGTraceIntegrator
         tr = mesh->GetFaceElementTransformations(iFace);
         IntegrationPointTransformation curTrans = tr->Elem1No == iCell ? tr->Loc1 : tr->Loc2;

         int intorder = (min(tr->Elem1->OrderW(), tr->Elem2->OrderW()) +
                        2*max(fes->GetFE(tr->Elem1No)->GetOrder(), fes->GetFE(tr->Elem2No)->GetOrder()));

         const IntegrationRule *ir = &IntRules.Get(tr->FaceGeom, intorder);
         // if (iCell == 1789 && myRank == 17) 
         // {
         //    std::cout << "tr->FaceGeom = " << tr->FaceGeom << ", intorder = " << intorder << std::endl;
         //    std::cout << "npoints = " << ir->GetNPoints() << std::endl;
         //    cout << "numFace = " << iFace << ";\n";
         // }

         updateYMin(*ir, &curTrans, elfun1_mat, iCell, yMin);

      } // for iFace

      // run through faces again to compute values in gauss points
      for (int iFace : shared_face_numbers)
      {
         // Integration order calculation from DGTraceIntegrator
         tr = mesh->GetSharedFaceTransformations(iFace);
         IntegrationPointTransformation curTrans = tr->Loc1;

         int intorder;
            intorder = tr->Elem1->OrderW() + 2*fes->GetFE(tr->Elem1No)->GetOrder();

         const IntegrationRule *ir = &IntRules.Get(tr->FaceGeom, intorder);

         // if (iCell == 1789 && myRank == 17) 
         // {
         //    std::cout << "tr->FaceGeom = " << tr->FaceGeom << ", intorder = " << intorder << std::endl;
         //    std::cout << "npoints = " << ir->GetNPoints() << std::endl;
         //    cout << "numFace = " << iFace << ";\n";
         // }
         updateYMin(*ir, &curTrans, elfun1_mat, iCell, yMin);
      } // for iFace

      // update ymin for vertices too
      // const IntegrationRule *irVertices = Geometries.GetVertices(fe->GetGeomType());

      // updateYMin(*irVertices, NULL, elfun1_mat, iCell, yMin);

      // if (iCell == 191)
      // {
      //    cout << "yMin = ";
      //    yMin.Print(cout);
      // }

      // for (int iPoint = 0; iPoint < irVertices->GetNPoints(); ++iPoint)
      // {
      //    const IntegrationPoint &ip = irVertices->IntPoint(iPoint);

      //    // Calculate basis functions on elements at the face
      //    fe->CalcShape(ip, shape);

      //    Vector funval1Vert(num_equation);

      //    // Interpolate elfun at the point
      //    elfun1_mat.MultTranspose(shape, funval1Vert);

      //    if (!StateIsPhysical(funval1Vert,dim))
      //    {
      //       for (int iEq = 0; iEq < num_equation; ++iEq)
      //          for (int iDof = 0; iDof < nDofs; ++iDof)  
      //             elfun1_mat(iDof, iEq) = uMean(iEq,iCell);
      //    }
      // }
      
      // find gauss points and loop through them
      //const IntegrationRule *ir = &IntRules.Get(el_trans->ElementType, fe->GetOrder());

      // replace values in gauss points for cell to the new one
      for (int iEq = 0; iEq < num_equation; ++iEq)
      {
         if (yMin[iEq] < 1.0)
           linearize(iCell, uMean, elfun1_mat);

         for (int iDof = 0; iDof < nDofs; ++iDof)
         {
            

            elfun1_mat(iDof, iEq) = uMean(iEq,iCell) + yMin[iEq] * (elfun1_mat(iDof, iEq) - uMean(iEq,iCell));
         }
      }

      // Last hope limiter
      const IntegrationRule *irVertices = Geometries.GetVertices(fe->GetGeomType());

      for (int iPoint = 0; iPoint < irVertices->GetNPoints(); ++iPoint)
      {
         const IntegrationPoint &ip = irVertices->IntPoint(iPoint);

         // Calculate basis functions on elements at the face
         fe->CalcShape(ip, shape);

         Vector funval1Vert(num_equation);

         // Interpolate elfun at the point
         elfun1_mat.MultTranspose(shape, funval1Vert);

         // double alpha = 1.0;
         // const double minEps = 1e-10;

         // if (!StateIsPhysical(funval1Vert,dim))
         // {
         //    // let's try to correct slope aсcording to positive internal energy
            
         //    alpha = min( alpha, 2.0 * funval1Vert[0] * (funval1Vert[dim+1] - minEps) / (funval1Vert[1]*funval1Vert[1] + funval1Vert[2]*funval1Vert[2]) );
         //    // for (int iEq = 0; iEq < num_equation; ++iEq)
         //    //    for (int iDof = 0; iDof < nDofs; ++iDof)  
         //    //       elfun1_mat(iDof, iEq) = uMean(iEq,iCell);
         // }

         // if (alpha < 1)
         // {
         //    // cout << "find alpha = " << alpha << " in cell #" << iCell << endl;
         //    for (int iEq = 1; iEq < num_equation - 1; ++iEq)
         //       for (int iDof = 0; iDof < nDofs; ++iDof) 
         //       {
         //          // const double beta = (sqrt(alpha) - 1.0) * uMean(iEq,iCell) / (elfun1_mat(iDof, iEq) - uMean(iEq,iCell)) + sqrt(alpha);
         //          // cout << " \tfind beta = " << beta << " for eqn #" << iEq << endl;
 
         //          // elfun1_mat(iDof, iEq) = uMean(iEq,iCell) + beta * (elfun1_mat(iDof, iEq) - uMean(iEq,iCell));
         //          elfun1_mat(iDof, iEq) = alpha * elfun1_mat(iDof, iEq);
         //       }
         // }
         if (!StateIsPhysical(funval1Vert,dim))
         {
            for (int iEq = 0; iEq < num_equation; ++iEq)
               for (int iDof = 0; iDof < nDofs; ++iDof)  
                  elfun1_mat(iDof, iEq) = uMean(iEq,iCell);
         }
      }

      // save limited solution values to the other vector
      xNew.SetSubVector(vdofs, el_x);
      
       // clean stencil
      internal_face_numbers.DeleteAll();
      shared_face_numbers.DeleteAll();
      stencil_num.DeleteAll();

   } // for iCell

   // replace solution values to the new one
   x = xNew;
   mesh->ExchangeFaceNbrData();
   parGridX.ExchangeFaceNbrData();

   // cout << "limited solution";
   // x.Print(cout);

   // cout << "end BJ" << endl;

};
