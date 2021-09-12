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
      Me_inv(vfes.GetNE() > 0 ? vfes.GetFE(0)->GetDof() : 0, vfes.GetNE() > 0 ? vfes.GetFE(0)->GetDof() : 0, vfes.GetNE()),
      state(num_equation),
      f(num_equation, dim),
      flux(vfes.GetNDofs(), dim, num_equation),
      z(A.Height())
{
    // Standard local assembly and inversion for energy mass matrices.
    int dof = 0;
    if (vfes.GetNE() > 0)
        dof = vfes.GetFE(0)->GetDof();
    DenseMatrix Me(dof);
    DenseMatrixInverse inv(&Me);
    MassIntegrator mi;
    for (int i = 0; i < vfes.GetNE(); i++)
    {
        mi.AssembleElementMatrix(*vfes.GetFE(i), *vfes.GetElementTransformation(i), Me);
        inv.Factor();
        inv.GetInverseMatrix(Me_inv(i));
        //std::cout << "ME_INV_" << i << std::endl;
        //Me_inv(i).Print(std::cout);
    }
}

void FE_Evolution::UpdateInverseMassMatrix()
{
    if (vfes.GetNE() > 0)
    {
        Me_inv.SetSize(vfes.GetFE(0)->GetDof(), vfes.GetFE(0)->GetDof(), vfes.GetNE());

        const int dof = vfes.GetFE(0)->GetDof();
        DenseMatrix Me(dof);
        DenseMatrixInverse inv(&Me);
        MassIntegrator mi;
        for (int i = 0; i < vfes.GetNE(); i++)
        {
            mi.AssembleElementMatrix(*vfes.GetFE(i), *vfes.GetElementTransformation(i), Me);
            inv.Factor();
            inv.GetInverseMatrix(Me_inv(i));
            //std::cout << "ME_INV_" << i << std::endl;
            //Me_inv(i).Print(std::cout);
        }
    }
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
    // 0. Reset wavespeed computation before operator application.
    max_char_speed = 0.;
  

    // cout << "IN FE_Evolution::Mult X = " << endl;
    // x.Print(cout);
    z.SetSize(x.Size());


    // 1. Create the vector z with the face terms -<F.n(u), [w]>.
    A.Mult(x, z);
// 
    // if (myRank == 6)
    // {
    // cout << "IN FE_Evolution::Mult z = " << endl;
    // z.Print(cout);
    // }

    // 2. Add the element terms.
    // i.  computing the flux approximately as a grid function by interpolating
    //      at the solution nodes.
    // ii. multiplying this grid function by a (constant) mixed bilinear form for
    //      each of the num_equation, computing (F(u), grad(w)) for each equation.

    DenseMatrix xmat(x.GetData(), vfes.GetNDofs(), num_equation);
    // cout << "IN FE_Evolution after dense mat " << endl;

    // cout << vfes.GetNDofs() << ' ' << dim << ' ' << num_equation << endl;

    // cout << flux.SizeI() << ' ' << flux.SizeJ() << ' ' << flux.SizeK() << endl;

    flux.SetSize(vfes.GetNDofs(), dim, num_equation);
     // cout << "IN FE_Evolution before get flux = " << endl;

    GetFlux(xmat, flux);
    // cout << "IN FE_Evolution after get flux = " << endl;
    // if (myRank == 6)
    // {
    // cout << " z before AddMult = " << endl;
    // z.Print(cout);
    // }

    for (int k = 0; k < num_equation; k++)
    {
        Vector fk(flux(k).GetData(), dim * vfes.GetNDofs());
        Vector zk(z.GetData() + k * vfes.GetNDofs(), vfes.GetNDofs());
        // if (myRank == 6)
        // {
        //     cout << " flux(k) before AddMult = " << endl;
        //     flux(k).Print(cout);
        //     cout << " fk before AddMult = " << endl;
        //     fk.Print(cout);
        //     cout << " zk before AddMult = " << endl;
        //     zk.Print(cout);
        // }
        Aflux->AddMult(fk, zk);
    }

    // 3. Multiply element-wise by the inverse mass matrices.
    Vector zval;
    Array<int> vdofs;
    int dof = 0;//vfes.GetNE() > 0 ? vfes.GetFE(0)->GetDof() : 0;
    DenseMatrix zmat, ymat(dof, num_equation);

    for (int i = 0; i < vfes.GetNE(); i++)
    {
        // Return the vdofs ordered byNODES
        vfes.GetElementVDofs(i, vdofs);
        z.GetSubVector(vdofs, zval);
        dof = vfes.GetFE(i)->GetDof();

        zmat.SetSize(dof, num_equation);
        ymat.SetSize(dof, num_equation);

        zmat.UseExternalData(zval.GetData(), dof, num_equation);
        // if (myRank == 6)
        // {
        //     std::cout << "zval" << endl;
        //     zval.Print(std::cout);
        // }
        mfem::Mult(Me_inv(i), zmat, ymat);
        
        y.SetSubVector(vdofs, ymat.GetData());

        // if (myRank == 6)
        // {
        //     std::cout << "zmat[" << i << "]  = ";
        //     zmat.Print(std::cout);
        // }
    }
    // cout << "IN FE_Evolution::end " << endl;

    // if (myRank == 6)
    // {
    // std::cout << "=======\n y after  = ";
    //     y.Print(std::cout);
    // }
}

// Compute the flux at solution nodes.
void FE_Evolution::GetFlux(const DenseMatrix &x, DenseTensor &flux) const
{
    const int dof = flux.SizeI();
    const int dim = flux.SizeJ();

    // if (myRank == 6) cout << "in GetFlux dof = " << dof << "; dim = " << dim << endl;

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

        // Update max char speed
        const double mcs = ComputeMaxCharSpeed(state, dim);
        if (mcs > max_char_speed) { max_char_speed = mcs; }
    }
}