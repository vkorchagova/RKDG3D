#include "fe_evolution.hpp"
#include "physics.hpp"


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(FiniteElementSpace &_vfes,
                              Operator &_A, SparseMatrix *_Aflux)
   : TimeDependentOperator(_A.Height()),
      dim(num_equation - 2),
      vfes(_vfes),
      A(_A),
      Aflux(_Aflux),
      Me_inv(vfes.GetNE()),
      state(num_equation),
      f(num_equation, dim),
      flux(vfes.GetNDofs(), dim, num_equation),
      z(A.Height())
{
   // Standard local assembly and inversion for energy mass matrices.
   int dof = 0;

   DenseMatrix Me(dof);
   DenseMatrixInverse inv(&Me);
   MassIntegrator mi;

   // double* vertex = new double[2];
   // for (int i = 0; i < vfes.GetNE(); i++)
   // {
   //    if (i == 50 || i == 1562)
   //    {   
   //    Array<int> indices;
   //    vfes.GetMesh()->GetElementVertices(i, indices);
   //    for (int iii = 0; iii < indices.Size(); iii++)
   //    {
   //         vertex = vfes.GetMesh()->GetVertex(indices[iii]);
   //         std::cout << std::std::setprecision(30) << "vertex[" << i << "][" << indices[iii] << "] : " << vertex[0] << ' ' << vertex[1] << std::endl;
   //    }
   // }
   // }
   
   for (int i = 0; i < vfes.GetNE(); i++)
   {
       // if (i == 50 || i == 1562)
       // { 
       //    Array<int> indices;
       //    vfes.GetMesh()->GetElementVertices(i, indices);
       //    for (int iii = 0; iii < indices.Size(); iii++)
       //    {
       //         vertex = vfes.GetMesh()->GetVertex(indices[iii]);
       //         std::cout << "vertex[" << i << "][" << indices[iii] << "] : " << vertex[0] << ' ' << vertex[1] << std::endl;
       //    }
       // }

       dof = vfes.GetFE(i)->GetDof();

       Me.SetSize(dof);
       mi.AssembleElementMatrix(*vfes.GetFE(i), *vfes.GetElementTransformation(i), Me);
       
       
       inv.Factor(Me);
       // inv.TestInversion();
       Me_inv[i].SetSize(dof,dof);
       inv.GetInverseMatrix(Me_inv[i]);

       // for (int ii = 0; ii < dof; ++ii)
       //    for (int jj = 0; jj < dof; ++jj)
       //    {
       //         if (fabs(Me_inv[i].Elem(ii,jj)) < 1e-16)
       //              Me_inv[i].Elem(ii,jj) = 0.0;
       //    }
       // if (i == 50 || i == 1562)
       // {
       //    Me_inv[i].Print(std::cout  << std::std::setprecision(30) << "Constr ME_INV_" << i << std::endl);
       // }
   }

   // for (int i = 50; i <= 1562; i+= 1512)
   // {   
   //    Array<int> indices;
   //    vfes.GetMesh()->GetElementVertices(i, indices);
   //    for (int iii = 0; iii < indices.Size(); iii++)
   //    {
   //         vertex = vfes.GetMesh()->GetVertex(indices[iii]);
   //         std::cout << "vertex[" << i << "][" << indices[iii] << "] : " << vertex[0] << ' ' << vertex[1] << std::endl;
   //    }
   // }
}

void FE_Evolution::UpdateInverseMassMatrix()
{
   if (vfes.GetNE() > 0)
   {
       Me_inv.resize(vfes.GetNE());
       const int dof = vfes.GetFE(0)->GetDof();
       DenseMatrix Me(dof);
       DenseMatrixInverse inv(&Me);
       MassIntegrator mi;
       for (int i = 0; i < vfes.GetNE(); i++)
       {
            mi.AssembleElementMatrix(*vfes.GetFE(i), *vfes.GetElementTransformation(i), Me);
            inv.Factor();
            Me_inv[i].SetSize(dof,dof);
            inv.GetInverseMatrix(Me_inv[i]);
       }
   }
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // 0. Reset wavespeed computation before operator application.
   max_char_speed = 0.;
  

   // std::cout << "IN FE_Evolution::Mult X = " << std::endl;
   // x.Print(std::cout);
   z.SetSize(x.Size());


   // 1. Create the vector z with the face terms -<F.n(u), [w]>.
   A.Mult(x, z);
// 
   // if (myRank == 6)
   // {
   // std::cout << "IN FE_Evolution::Mult z = " << std::endl;
   // z.Print(std::cout);
   // }
   Vector zval;
   Array<int> vdofs;
   // for (int i = 0; i < vfes.GetNE(); i++)
   // {
   //    vfes.GetElementVDofs(i, vdofs);
   //    z.GetSubVector(vdofs, zval);
   //    if (i == 0 || i == 1512)
   //    { 
   //         zval.Print(std::cout << std::std::setprecision(18) << "zval[" << i << "] after fluxes = " << std::endl);
   //    }
   // }

   // 2. Add the element terms.
   // i.  computing the flux approximately as a grid function by interpolating
   //          at the solution nodes.
   // ii. multiplying this grid function by a (constant) mixed bilinear form for
   //      each of the num_equation, computing (F(u), grad(w)) for each equation.

   DenseMatrix xmat(x.GetData(), vfes.GetNDofs(), num_equation);
   // std::cout << "IN FE_Evolution after dense mat " << std::endl;

   // std::cout << vfes.GetNDofs() << ' ' << dim << ' ' << num_equation << std::endl;

   // std::cout << flux.SizeI() << ' ' << flux.SizeJ() << ' ' << flux.SizeK() << std::endl;

   flux.SetSize(vfes.GetNDofs(), dim, num_equation);
    // std::cout << "IN FE_Evolution before get flux = " << std::endl;

   GetFlux(xmat, flux);
   // std::cout << "IN FE_Evolution after get flux = " << std::endl;
   // if (myRank == 6)
   // {
   // std::cout << " z before AddMult = " << std::endl;
   // z.Print(std::cout);
   // }

   for (int k = 0; k < num_equation; k++)
   {
       Vector fk(flux(k).GetData(), dim * vfes.GetNDofs());
       Vector zk(z.GetData() + k * vfes.GetNDofs(), vfes.GetNDofs());
       // if (myRank == 6)
       // {
            // std::cout << " flux(k) before AddMult = " << std::endl;
            // flux(k).Print(std::cout);
            // std::cout << " fk before AddMult = " << std::endl;
            // fk.Print(std::cout);
            // std::cout << " zk before AddMult = " << std::endl;
            // zk.Print(std::cout);
       // }
       Aflux->AddMult(fk, zk);
   }

   // 3. Multiply element-wise by the inverse mass matrices.
   
   int dof = 0;//vfes.GetNE() > 0 ? vfes.GetFE(0)->GetDof() : 0;
   DenseMatrix zmat, ymat(dof, num_equation);

   for (int i = 0; i < vfes.GetNE(); i++)
   {
       // Return the vdofs ordered byNODES
       vfes.GetElementVDofs(i, vdofs);
       z.GetSubVector(vdofs, zval);
       dof = vfes.GetFE(i)->GetDof();

       // zmat.SetSize(dof, num_equation);
       ymat.SetSize(dof, num_equation);

       zmat.Reset(zval.GetData(), dof, num_equation);
       // if (i == 0 || i == 1512)
       // {
       //    zval.Print(std::cout << std::std::setprecision(18) << "zval[" << i << ']' << std::endl);
       // }
       mfem::Mult(Me_inv[i], zmat, ymat);

       
       // if (i == 50 || i == 1562)
       // {
       //    Me_inv[i].Print(std::cout  << std::std::setprecision(30) << "ME_INV_" << i << std::endl);
       // }
       
       y.SetSubVector(vdofs, ymat.GetData());

       // if (i == 0 || i == 1512)
       // {
       //    zmat.Print(std::cout << std::std::setprecision(18) << "zmat[" << i << "]  = ");
       //    ymat.Print(std::cout << std::std::setprecision(18) << "ymat[" << i << "]  = ");
       // }
   }
   // std::cout << "IN FE_Evolution::end " << std::endl;

   // if (myRank == 6)
   // {
   // std::cout << "=======\n y after  = ";
   //    y.Print(std::cout);
   // }
}

// Compute the flux at solution nodes.
void FE_Evolution::GetFlux(const DenseMatrix &x, DenseTensor &flux) const
{
   const int dof = flux.SizeI();
   const int dim = flux.SizeJ();

   // if (myRank == 6) std::cout << "in GetFlux dof = " << dof << "; dim = " << dim << std::endl;

   for (int i = 0; i < dof; i++)
   {
       for (int k = 0; k < num_equation; k++) 
            state(k) = x(i, k);

       ComputeFlux(state, dim, f);

       for (int d = 0; d < dim; d++)
       {
            for (int k = 0; k < num_equation; k++)
            {
                flux(i, d, k) = f(k, d);
            }
       }

       // // Update max char speed
       // const double mcs = ComputeMaxCharSpeed(state, dim);
       // if (mcs > max_char_speed) { max_char_speed = mcs; }
   }
}