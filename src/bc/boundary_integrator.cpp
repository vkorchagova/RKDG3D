#include "boundary_integrator.hpp"

// Implementation of class BoundaryIntegrator
BoundaryIntegrator::BoundaryIntegrator(RiemannSolver &rsolver_, const int dim) :
   rsolver(rsolver_),
   shape1(num_equation),
   funval1(num_equation),
   funval2(num_equation),
   nor(dim),
   fluxN(num_equation),
   dim(dim) {}

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

   // Vector eip11(3);

   DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equation);

   DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equation);

   // std::cout << "======= elfun\n";
   // elfun.Print(std::cout);
   // if ((Tr.Elem1No == 0) && myRank == 3) 
   // {
   //   std::cout << "======= elfun1_mat\n";
   //   elfun1_mat.Print(std::cout);
   //   std::cout << "======= elfun2_mat\n";
   //   elfun2_mat.Print(std::cout);
   // }

   // Integration order calculation from DGTraceIntegrator
   int intorder;
   if (Tr.Elem2No >= 0)
      intorder = (std::min(Tr.Elem1->OrderW(), Tr.Elem2->OrderW()) +
                   2*std::max(el1.GetOrder(), el2.GetOrder()));
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
   //   IntegrationRule nds1 = el1.GetNodes();
   //   IntegrationRule nds2 = el2.GetNodes();

   //   for (int i = 0; i < nds1.GetNPoints(); i++)
   //   {
   //       const IntegrationPoint &ip = nds1.IntPoint(i);
   //       std::cout << "node cell1 = " << ip.x << " " << ip.y << std::endl;
   //   }
   //   for (int i = 0; i < nds2.GetNPoints(); i++)
   //   {
   //       const IntegrationPoint &ip = nds2.IntPoint(i);
   //       std::cout << "node cell2 = " << ip.x << " " << ip.y << std::endl;
   //   }
   // }

   //std::cout << "Tr.FaceGeom = " << Tr.FaceGeom << ", intorder = " << intorder << std::endl;
   // std::cout << "npoints = " << ir->GetNPoints() << std::endl;


   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      // std::cout << "ipoint = " << ip.index << " " << ip.x << " " << ip.y << " " << ip.z << std::endl;

      // Tr.Loc1.Transform(ip, eip1);
      // /// ========= 
      // // shape.SetSize(FElem->GetDof());
      // // trans.SetSize(PointMat.Height());
      // Tr.Loc1.Transf.GetFE()->CalcShape(ip, shape1);
      // DenseMatrix PointMat = Tr.Loc1.Transf.GetPointMat();
      // PointMat.Mult(shape1, eip11);
      Tr.SetAllIntPoints(&ip); // set face and element int. points

      // Calculate basis functions on both elements at the face
      el1.CalcShape(Tr.GetElement1IntPoint(), shape1);

      // for (int ii = 0; ii < shape1.Size(); ++ii)
      // {
      //   if (fabs(shape1(ii)) < 1e-16)
      //       shape1(ii) = 0.0;
      // }

      // if (Tr.Elem1No == 0 || Tr.Elem1No == 1512)
      // {
      //   std::cout << "Get RS info for Tr.Elem1No = " << Tr.Elem1No << std::endl;
      //   // std::cout << "Tr.FaceGeom = " << Tr.FaceGeom << std::endl;
      //   std::cout  << std::std::setprecision(20) << "ipoint = " << ip.index << " " << ip.x << " " << ip.y << " " << ip.z << std::endl;
      //   // std::cout  << std::std::setprecision(20) << "gp_eip1 = " << eip1.x << " " << eip1.y <<  std::endl;
      //   // Tr.Loc1.Transf.GetPointMat().Print(std::cout << "PointMat = " << std::std::setprecision(20));
      //   // eip11.Print(std::cout  << std::std::setprecision(20) << "gp_eip11 = " );
      //   shape1.Print(std::cout << "shape1 = " << std::std::setprecision(20));
      //   // elfun1_mat.Print(std::cout << "elfun1_mat = " << std::std::setprecision(20));
         
      // }

      // Get the normal vector and the flux on the face
      CalcOrtho(Tr.Face->Jacobian(), nor);
    
      double normag = 0;
      for (int i = 0; i < nor.Size(); i++)
      {
         normag += nor(i) * nor(i);
      }
      normag = sqrt(normag);
   
      nor *= 1.0/nor.Norml2();

      //std::cout << "el1.Space = " << el1.Space() << std::endl;
      //std::cout << "el2.Space = " << el2.Space() << std::endl;

      // std::cout << "el1.CalcShape = "; shape1.Print(std::cout);
      // std::cout << "el2.CalcShape = "; shape2.Print(std::cout);
 
      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(shape1, funval1);
      
      computeRightState(funval1, funval2, nor);

      // if ((Tr.Elem1No == 0 || Tr.Elem2No == 0) ) 
      // {
      //   std::cout << " 0 funval1 = "; funval1.Print(std::cout);
      //   std::cout << " 0 funval2 = "; funval2.Print(std::cout);
      // }

      // if ((Tr.Elem1No == 966 || Tr.Elem2No == 966) ) 
      // {
      //   std::cout << " 966 funval1 = "; funval1.Print(std::cout);
      //   std::cout << " 966 funval2 = "; funval2.Print(std::cout);
      // }

      // Tr.Face->SetIntPoint(&ip);

      bool debugRS = false; 

      // if (Tr.Elem1No == 0 || Tr.Elem1No == 1512)
      // {
      //   debugRS = true;
       //  }
 

      const double mcs = rsolver.Eval(funval1, funval2, nor, fluxN, debugRS);

      // if (Tr.Elem1No == 0 || Tr.Elem1No == 1512)
      //  {
      //   std::cout << "After RS Compute - " << Tr.Elem1No << std::endl;
      //   std::cout << "---" << std::endl;
      // }

      if (mcs < 0) 
      {
         std::cout << "Number of neighbours: " << Tr.Elem1No << ' ' << Tr.Elem2No << std::endl;
         exit(1);
      }

      // std::cout << "funval1: ";
      // funval1.Print(std::cout);
      // std::cout << "funval2: ";
      // funval2.Print(std::cout);
      // std::cout << "\tmcs = " << mcs << std::endl;
      // std::cout << "\tfluxN = ";
      // fluxN.Print(std::cout);


      // Update max char speed
      if (mcs > max_char_speed) { max_char_speed = mcs; }

      fluxN *= ip.weight;
      // std::cout << "nor = ";
      // nor.Print(std::cout);
      // std::cout << "ip weight = " << ip.weight << std::endl;
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
