#include "domain_integrator.hpp"

using namespace std;

// Implementation of class DomainIntegrator
DomainIntegrator::DomainIntegrator(const int dim) /*: flux(num_equation, dim) */{ }

void DomainIntegrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                                   const FiniteElement &test_fe,
                                                   ElementTransformation &Tr,
                                                   DenseMatrix &elmat)
{
   // Assemble the form (vec(v), grad(w))

   // Trial space = vector L2 space (mesh dim)
   // Test space  = scalar L2 space

   const int dof_trial = trial_fe.GetDof();
   const int dof_test = test_fe.GetDof();
   const int dim = trial_fe.GetDim();

   shape.SetSize(dof_trial);
   dshapedr.SetSize(dof_test, dim);
   dshapedx.SetSize(dof_test, dim);

   elmat.SetSize(dof_test, dof_trial * dim);
   elmat = 0.0;

   const int maxorder = std::max(trial_fe.GetOrder(), test_fe.GetOrder());
   const int intorder = 2 * maxorder;
   const IntegrationRule *ir = &IntRules.Get(trial_fe.GetGeomType(), intorder);

   for (int i = 0; i < ir->GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      // Calculate the shape functions
      trial_fe.CalcShape(ip, shape);
      shape *= ip.weight;

      // Compute the physical gradients of the test functions
      Tr.SetIntPoint(&ip);
      test_fe.CalcDShape(ip, dshapedr);
      Mult(dshapedr, Tr.AdjugateJacobian(), dshapedx);

      for (int d = 0; d < dim; d++)
      {
         for (int j = 0; j < dof_test; j++)
         {
            for (int k = 0; k < dof_trial; k++)
            {
                elmat(j, k + d * dof_trial) += shape(k) * dshapedx(j, d);
            }
         }
      }
   }
}


// void DomainIntegrator::ComputeElementFlux
// ( const FiniteElement &el, ElementTransformation &Trans,
//   Vector &u, const FiniteElement &fluxelem, Vector &flux, bool with_coef )
// {
//   int i, j, nd, dim, spaceDim, fnd;

//   nd = el.GetDof();
//   dim = el.GetDim();
//   spaceDim = Trans.GetSpaceDim();

//   if (VQ)
//   {
//       MFEM_VERIFY(VQ->GetVDim() == spaceDim,
//                     "Unexpected dimension for VectorCoefficient");
//   }
//   if (MQ)
//   {
//       MFEM_VERIFY(MQ->GetWidth() == spaceDim,
//                     "Unexpected width for MatrixCoefficient");
//       MFEM_VERIFY(MQ->GetHeight() == spaceDim,
//                     "Unexpected height for MatrixCoefficient");
//   }

// #ifdef MFEM_THREAD_SAFE
//   DenseMatrix dshape(nd,dim), invdfdx(dim, spaceDim);
//   DenseMatrix M(MQ ? spaceDim : 0);
//   Vector D(VQ ? VQ->GetVDim() : 0);
// #else
//   dshape.SetSize(nd,dim);
//   invdfdx.SetSize(dim, spaceDim);
//   M.SetSize(MQ ? spaceDim : 0);
//   D.SetSize(VQ ? VQ->GetVDim() : 0);
// #endif
//   vec.SetSize(dim);
//   vecdxt.SetSize(spaceDim);
//   pointflux.SetSize(MQ ? spaceDim : 0);

//   const IntegrationRule &ir = fluxelem.GetNodes();
//   fnd = ir.GetNPoints();
//   flux.SetSize( fnd * spaceDim );

//   for (i = 0; i < fnd; ++i)
//   {
//       const IntegrationPoint &ip = ir.IntPoint(i);
//       el.CalcDShape(ip, dshape);
//       dshape.MultTranspose(u, vec);

//       Trans.SetIntPoint (&ip);
//       CalcInverse(Trans.Jacobian(), invdfdx);
//       invdfdx.MultTranspose(vec, vecdxt);

//       if (!MQ && !VQ)
//       {
//          if (Q && with_coef)
//          {
//              vecdxt *= Q->Eval(Trans,ip);
//          }
//          for (j = 0; j < spaceDim; j++)
//          {
//              flux(fnd*j+i) = vecdxt(j);
//          }
//       }
//       else
//       {
//          if (MQ)
//          {
//              MQ->Eval(M, Trans, ip);
//              M.Mult(vecdxt, pointflux);
//          }
//          else
//          {
//              VQ->Eval(D, Trans, ip);
//              for (int j=0; j<spaceDim; ++j)
//              {
//                pointflux[j] = D[j] * vecdxt[j];
//              }

//          }
//          for (j = 0; j < spaceDim; j++)
//          {
//              flux(fnd*j+i) = pointflux(j);
//          }
//       }
//   }
// }

// double DomainIntegrator::ComputeFluxEnergy
// ( const FiniteElement &fluxelem, ElementTransformation &Trans,
//   Vector &flux, Vector* d_energy)
// {
//   int nd = fluxelem.GetDof();
//   int dim = fluxelem.GetDim();
//   int spaceDim = Trans.GetSpaceDim();

// #ifdef MFEM_THREAD_SAFE
//   DenseMatrix M;
// #endif

//   shape.SetSize(nd);
//   pointflux.SetSize(spaceDim);
//   if (d_energy) { vec.SetSize(spaceDim); }
//   if (MQ) { M.SetSize(spaceDim); }

//   int order = 2 * fluxelem.GetOrder(); // <--
//   const IntegrationRule *ir = &IntRules.Get(fluxelem.GetGeomType(), order);

//   double energy = 0.0;
//   if (d_energy) { *d_energy = 0.0; }

//   for (int i = 0; i < ir->GetNPoints(); ++i)
//   {
//       const IntegrationPoint &ip = ir->IntPoint(i);
//       fluxelem.CalcShape(ip, shape);

//       pointflux = 0.0;
//       for (int k = 0; k < spaceDim; k++)
//       {
//          for (int j = 0; j < nd; j++)
//          {
//              pointflux(k) += flux(k*nd+j)*shape(j);
//          }
//       }

//       Trans.SetIntPoint(&ip);
//       double w = Trans.Weight() * ip.weight;

//       if (!MQ)
//       {
//          double e = (pointflux * pointflux);
//          if (Q) { e *= Q->Eval(Trans, ip); }
//          energy += w * e;
//       }
//       else
//       {
//          MQ->Eval(M, Trans, ip);
//          energy += w * M.InnerProduct(pointflux, pointflux);
//       }

//       if (d_energy)
//       {
//          // transform pointflux to the ref. domain and integrate the components
//          Trans.Jacobian().MultTranspose(pointflux, vec);
//          for (int k = 0; k < dim; k++)
//          {
//              (*d_energy)[k] += w * vec[k] * vec[k];
//          }
//          // TODO: Q, MQ
//       }
//   }

//   return energy;
// }
