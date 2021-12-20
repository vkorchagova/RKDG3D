#include "check_total_energy.hpp"

double ComputeTotalEnergy(ParMesh* mesh, ParFiniteElementSpace* vfes, ParGridFunction& sol)
{
    double totalEnergy = 0.0;
    double myTotalEnergy = 0.0;
    double totalEnergy_el = 0.0;

    const FiniteElement *fe;
    Vector el_sol;
    Array<int> vdofs;
    ElementTransformation* el_trans;
    Vector nor;

    for (int iCell = 0; iCell < vfes->GetNE(); ++iCell)
    {
        totalEnergy_el = 0.0;
        fe = vfes->GetFE(iCell);
        const int nDofs = fe->GetDof();
        nor.SetSize(nDofs);

        vfes->GetElementVDofs(iCell, vdofs);
        sol.GetSubVector(vdofs, el_sol);

        el_trans = mesh->GetElementTransformation(iCell);


        // el_sol.Print(cout << "el_sol for cell #" << iCell << ": ");

        const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), fe->GetOrder() + 1);

        for (int iDof = 0; iDof < nDofs; ++iDof)
        {
            const IntegrationPoint ip = ir->IntPoint(iDof);
            // el_trans->SetIntPoint(&ip);
            // CalcOrtho(el_trans->Jacobian(), nor);
            

            totalEnergy_el += ip.weight * el_sol[(num_equation - 1)*nDofs + iDof];
        }

        // cout << el_trans->Jacobian().Det() << ' ' << totalEnergy_el << ' ' << totalEnergy_el * el_trans->Jacobian().Det() << endl;

        // el_trans->Jacobian().Print(cout << "Jacobian for cell #" << iCell << "= ");

        myTotalEnergy += totalEnergy_el * el_trans->Jacobian().Det();


    }

    MPI_Reduce(
                &myTotalEnergy,
                &totalEnergy,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                0,
                mesh->GetComm()
            );


    return totalEnergy;
}

// void ComputeU(ParGridFunction& sol, ParFiniteElementSpace* vfes, ParGridFunction& U, ParFiniteElementSpace* dfes)
// {
    
//     for (int i = 0; i < dfes.GetNE(); i++)
//     {
        
//     }
// }