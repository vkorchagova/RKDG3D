#include "limiter_multiplier.hpp"
#include "physics.hpp"

// bool pleaseWriteMe = false;

LimiterMultiplier::LimiterMultiplier(Indicator& _ind, Averager& _avgr, ParFiniteElementSpace*_fes, const Array<int>& _offsets, int _d) : Limiter(_ind, _avgr, _fes, _offsets, _d) 
{

};


void LimiterMultiplier::limit(const int iCell, const Vector& el_ind, DenseMatrix& elfun1_mat) 
{
   // get FE
   fe = fes->GetFE(iCell);

   const int nDofs = fe->GetDof();

   averager.readElementAverageByNumber(iCell, el_uMean);

   // replace solution to mean values
   for (int iEq = 0; iEq < num_equation; ++iEq)
      for (int iDof = 0; iDof < nDofs; ++iDof)
      {
         // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval before lim = " << elfun1_mat(iDof, iEq) << endl;}
         if (el_ind[iEq] < 1.0)
            linearize(iCell, el_uMean, elfun1_mat);

         for (int iDof = 0; iDof < nDofs; ++iDof)
         {
            elfun1_mat(iDof, iEq) = el_uMean[iEq] + el_ind[iEq] * (elfun1_mat(iDof, iEq) - el_uMean[iEq]);
         }
         // if (myRank == 0 && iCell == 0 && iEq == 0) {cout << "   funval after lim = " << elfun1_mat(iDof, iEq) << endl;}
      }

   // x.SetSubVector(el_vdofs, el_x);

   // // Last hope limiter
   // const IntegrationRule *irVertices = Geometries.GetVertices(fe->GetGeomType());

   // for (int iPoint = 0; iPoint < irVertices->GetNPoints(); ++iPoint)
   // {
   //    const IntegrationPoint &ip = irVertices->IntPoint(iPoint);

   //    // Calculate basis functions on elements at the face
   //    fe->CalcShape(ip, shape);

   //    Vector funval1Vert(num_equation);

   //    // Interpolate elfun at the point
   //    elfun1_mat.MultTranspose(shape, funval1Vert);

   //    // double alpha = 1.0;
   //    // const double minEps = 1e-10;

   //    // if (!StateIsPhysical(funval1Vert,dim))
   //    // {
   //    //    // let's try to correct slope aсcording to positive internal energy
         
   //    //    alpha = min( alpha, 2.0 * funval1Vert[0] * (funval1Vert[dim+1] - minEps) / (funval1Vert[1]*funval1Vert[1] + funval1Vert[2]*funval1Vert[2]) );
   //    //    // for (int iEq = 0; iEq < num_equation; ++iEq)
   //    //    //    for (int iDof = 0; iDof < nDofs; ++iDof)  
   //    //    //       elfun1_mat(iDof, iEq) = uMean(iEq,iCell);
   //    // }

   //    // if (alpha < 1)
   //    // {
   //    //    // cout << "find alpha = " << alpha << " in cell #" << iCell << endl;
   //    //    for (int iEq = 1; iEq < num_equation - 1; ++iEq)
   //    //       for (int iDof = 0; iDof < nDofs; ++iDof) 
   //    //       {
   //    //          // const double beta = (sqrt(alpha) - 1.0) * uMean(iEq,iCell) / (elfun1_mat(iDof, iEq) - uMean(iEq,iCell)) + sqrt(alpha);
   //    //          // cout << " \tfind beta = " << beta << " for eqn #" << iEq << endl;

   //    //          // elfun1_mat(iDof, iEq) = uMean(iEq,iCell) + beta * (elfun1_mat(iDof, iEq) - uMean(iEq,iCell));
   //    //          elfun1_mat(iDof, iEq) = alpha * elfun1_mat(iDof, iEq);
   //    //       }
   //    // }
   //    if (!StateIsPhysical(funval1Vert,dim))
   //    {
   //       for (int iEq = 0; iEq < num_equation; ++iEq)
   //          for (int iDof = 0; iDof < nDofs; ++iDof)  
   //             elfun1_mat(iDof, iEq) = uMean(iEq,iCell);
   //    }
   // }

   // save limited solution values to the other vector
  // xNew.SetSubVector(el_vdofs, el_x);
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

   // std::cout << "elfun1_mat = ";
   // elfun1_mat.Print(std::cout);

   // get center of finite element
   mesh->GetElementCenter(iCell, el_center_phys);
   // if (iCell == 6950) el_center_phys.Print(std::cout);

   // move it to the element space
   el_trans->TransformBack(el_center_phys, el_center_ref);
   // if (iCell == 6950) std::cout << "el_center_ref = " << el_center_ref.index << " " << el_center_ref.x << " " << el_center_ref.y << " " << el_center_ref.z << std::endl;

   // compute gradient of el_shape functions in the center (OK for p = 1 because of 1,x,y,xy)
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

         elfun1_mat(iDof, iEq) = el_uMean(iEq) - \
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
