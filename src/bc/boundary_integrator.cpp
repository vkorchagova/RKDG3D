#include "boundary_integrator.hpp"

// Implementation of class BoundaryIntegrator
BoundaryIntegrator::BoundaryIntegrator(RiemannSolver &rsolver_, const int dim) :
   rsolver(rsolver_),
   funval1(num_equation),
   funval2(num_equation),
   nor(dim),
   fluxN(num_equation),
   dim(dim) { }

void BoundaryIntegrator::AssembleFaceVector(const FiniteElement &el1,
                                        const FiniteElement &el2,
                                        FaceElementTransformations &Tr,
                                        const Vector &elfun, Vector &elvect)
{
   // Compute the term <F.n(u),[w]> on the interior faces.
   const int dof1 = el1.GetDof();
   shape1.SetSize(dof1);

   elvect.SetSize((dof1) * num_equation);
   elvect = 0.0;

   DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equation);

   DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equation);

   // std::cout << "======= elfun\n";
   // elfun.Print(std::cout);
   // if ((Tr.Elem1No == 0) && myRank == 3) 
   // {
   //    std::cout << "======= elfun1_mat\n";
   //    elfun1_mat.Print(std::cout);
   //    std::cout << "======= elfun2_mat\n";
   //    elfun2_mat.Print(std::cout);
   // }

   // Integration order calculation from DGTraceIntegrator
   int intorder;
   if (Tr.Elem2No >= 0)
      intorder = (min(Tr.Elem1->OrderW(), Tr.Elem2->OrderW()) +
                  2*max(el1.GetOrder(), el2.GetOrder()));
   else
   {
      intorder = Tr.Elem1->OrderW() + 2*el1.GetOrder();
      // Tr.Elem2No = Tr.Elem1No;
      // Tr.Loc2 = Tr.Loc1;
   }
   if (el1.Space() == FunctionSpace::Pk)
   {
      intorder++;
   }
   //std::cout << "Tr.Elem1->OrderW() = " << Tr.Elem1->OrderW() << ", el1.GetOrder() = " << el1.GetOrder() << std::endl;

   const IntegrationRule *ir = &IntRules.Get(Tr.FaceGeom, intorder);

   // if ((Tr.Elem1No == 0) && myRank == 3) 
   // {
   //    IntegrationRule nds1 = el1.GetNodes();
   //    IntegrationRule nds2 = el2.GetNodes();

   //    for (int i = 0; i < nds1.GetNPoints(); i++)
   //    {
   //       const IntegrationPoint &ip = nds1.IntPoint(i);
   //       cout << "node cell1 = " << ip.x << " " << ip.y << endl;
   //    }
   //    for (int i = 0; i < nds2.GetNPoints(); i++)
   //    {
   //       const IntegrationPoint &ip = nds2.IntPoint(i);
   //       cout << "node cell2 = " << ip.x << " " << ip.y << endl;
   //    }
   // }

   //std::cout << "Tr.FaceGeom = " << Tr.FaceGeom << ", intorder = " << intorder << std::endl;
   // std::cout << "npoints = " << ir->GetNPoints() << std::endl;


   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      // std::cout << "ipoint = " << ip.index << " " << ip.x << " " << ip.y << " " << ip.z << std::endl;

      Tr.Loc1.Transform(ip, eip1);

      // Get the normal vector and the flux on the face
      CalcOrtho(Tr.Face->Jacobian(), nor);
     
      double normag = 0;
      for (int i = 0; i < nor.Size(); i++)
      {
         normag += nor(i) * nor(i);
      }
      normag = sqrt(normag);

      nor *= 1.0/nor.Norml2();

      // if ((Tr.Elem1No == 0 || Tr.Elem2No == 0) && myRank == 3) 
      // {
      //    Vector gp_phys1(2);
      //    Vector gp_phys2(2);
         
      //    // Tr.GetElement1Transformation().Transform(eip1,gp_phys1);
      //    // Tr.GetElement2Transformation().Transform(eip2,gp_phys2);
      //    std::cout << "face between elems: " << Tr.Elem1No << "|" << Tr.Elem2No << ";\n";
      //    std::cout << "\tgp_eip1 = " << eip1.x<< " " << eip1.y<<  endl;
      //    std::cout << "\tgp_eip2 = " << eip2.x << " " << eip2.y <<  endl;
      //    // std::cout << "\tgp_eltrans1 = " << gp_phys1[0] << " " << gp_phys1[1] <<  endl;
      //    // std::cout << "\tgp_eltrans2 = " << gp_phys2[0] << " " << gp_phys2[1] <<  endl;
      // }




      // if (myRank == 3 &&  Tr.Elem1No == 0) 
      // {
      //    Vector el_gp_phys(dim);
      //    el_trans = mesh->GetElementTransformation(iCell);
      //    el_trans->Transform(eip1,el_gp_phys);

      //    std::cout << "iFace = " << iFace << "; gpoint phys = " << el_gp_phys[0] << " " << el_gp_phys[1]   << std::endl;
      // }

      //Tr.Loc1.Transf.GetPointMat().Print(std::cout);
      //Tr.Loc2.Transf.GetPointMat().Print(std::cout);

      //std::cout << "e1point = " << eip1.index << " " << eip1.x << " " << eip1.y << " " << eip1.z << std::endl;
      //std::cout << "e2point = " << eip2.index << " " << eip2.x << " " << eip2.y << " " << eip2.z << std::endl;

      // Calculate basis functions on both elements at the face
      el1.CalcShape(eip1, shape1);

      //std::cout << "el1.Space = " << el1.Space() << std::endl;
      //std::cout << "el2.Space = " << el2.Space() << std::endl;

      // std::cout << "el1.CalcShape = "; shape1.Print(std::cout);
      // std::cout << "el2.CalcShape = "; shape2.Print(std::cout);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(shape1, funval1);
      
      computeRightState(funval1, funval2, nor);

      // if ((Tr.Elem1No == 0 || Tr.Elem2No == 0) && myRank == 3) 
      // {
      //    std::cout << "funval1 = "; funval1.Print(std::cout);
      //    std::cout << "funval2 = "; funval2.Print(std::cout);
      // }

      Tr.Face->SetIntPoint(&ip);

      

      const double mcs = rsolver.Eval(funval1, funval2, nor, fluxN);

      // cout << "funval1: ";
      // funval1.Print(cout);
      // cout << "funval2: ";
      // funval2.Print(cout);
      // cout << "\tmcs = " << mcs << endl;
      // std::cout << "\tfluxN = ";
      // fluxN.Print(std::cout);


      // Update max char speed
      if (mcs > max_char_speed) { max_char_speed = mcs; }

      fluxN *= ip.weight;
      // cout << "nor = ";
      // nor.Print(cout);
      // cout << "ip weight = " << ip.weight << endl;
      for (int k = 0; k < num_equation; k++)
      {
         for (int s = 0; s < dof1; s++)
         {
            elvect1_mat(s, k) -= fluxN(k) * shape1(s) * normag;
         }
      }
   }

   // std::cout << "elvect in the end\n";
   // elvect.Print(std::cout);
}
