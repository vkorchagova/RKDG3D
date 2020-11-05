#include "rk_explicit_limited.hpp"


ExplicitRKLimitedSolver::ExplicitRKLimitedSolver(int _s, const double *_a, const double *_b,
                                   const double *_c, Limiter& _l) : 
   s (_s),
   a (_a),
   b (_b),
   c (_c),
   limiter(_l)
{
   
   k = new Vector[s];
}

void ExplicitRKLimitedSolver::Init(TimeDependentOperator &_f)
{
   ODESolver::Init(_f);
   int n = f->Width();
   y.SetSize(n, mem_type);
   for (int i = 0; i < y.Size(); ++i) y(i) = 0.0;
   //if (myRank == 0) y.Print(cout);
   for (int i = 0; i < s; i++)
   {
      k[i].SetSize(n, mem_type);
   }
}


void ExplicitRKLimitedSolver::Step(Vector &x, double &t, double &dt)
{
   // std::cout << "=======\n X before limiting = ";
   //    x.Print(std::cout);
   // for (int iii = 0; iii < x.Size(); ++iii)
   //    if (x[iii] != x[iii])
   //    {
   //       cout << "Find NaN for x[" << iii << "], processor " << myRank << endl;
   //       //cout << "; iCell = " << iCell << "; tr elem 1 no = "<< tr->Elem1No << "; tr elem 2 no =" << tr->Elem2No << endl;
   //    }

   f->SetTime(t);
   f->Mult(x, k[0]);
   // for (int iii = 0; iii < k[0].Size(); ++iii)
   //    if (k[0][iii] != k[0][iii])
   //    {
   //       cout << "Find NaN for k[0][" << iii << "], processor " << myRank << endl;
   //       //cout << "; iCell = " << iCell << "; tr elem 1 no = "<< tr->Elem1No << "; tr elem 2 no =" << tr->Elem2No << endl;
   //    }
   // std::cout << "=======\n k[0]  = ";
   //    k[0].Print(std::cout);
   for (int l = 0, i = 1; i < s; i++)
   {
      //**/cout << "1/l|i|s = " << l << '|' << i << '|' << s << endl;
      add(x, a[l++]*dt, k[0], y);
      //**/cout << "2/l|i|s = " << l << '|' << i << '|' << s << endl;
      for (int j = 1; j  < i; j++)
      {
         y.Add(a[l++]*dt, k[j])                                                                                                                                                                                                                                                ;
         //**/cout << "3/l|i|s|j = " << l << '|' << i << '|' << s << '|' << j << endl;
      }

      // if (myRank == 3 || myRank == 0)
      // {
      //    std::cout << "Y after 1st stage before limit | ";
      //    std::cout << "rank = " << myRank << ", x = " << x[0] << std::endl;
      // }

      
      // limit y and compute rhs with good y
       limiter.limit(y);
      /// HERE WE NEED MPI& ???
      // std::cout << "=======\n X after 1st stage = ";
      // y.Print(std::cout);
      // std::cout << "next stage" << endl;
      // compute next rhs
      //  if (myRank == 3 || myRank == 0)
      // {
      //    std::cout << "Y after 1st stage after limit | ";
      //    std::cout << "rank = " << myRank << ", x = " << x[0] << std::endl;
      // }

      f->SetTime(t + c[i-1]*dt);
      f->Mult(y, k[i]);

      // for (int iii = 0; iii < k[i].Size(); ++iii)
      // if (k[i][iii] != k[i][iii])
      // {
      //    cout << "Find NaN for k[1][" << iii << "], processor " << myRank << endl;
      //    //cout << "; iCell = " << iCell << "; tr elem 1 no = "<< tr->Elem1No << "; tr elem 2 no =" << tr->Elem2No << endl;
      // }

      // std::cout << "=======\n k[i]  = ";
      // k[i].Print(std::cout);
   }

   // if (myRank == 3 || myRank == 0)
   // {
   //    std::cout << "X before 2nd stage before limit = ";
   //    std::cout << "rank = " << myRank << ", x = " << x[0] << std::endl;
   // }

   for (int i = 0; i < s; i++)
   {
      x.Add(b[i]*dt, k[i]);
   }

   // if (myRank == 3 || myRank == 0)
   // {
   //    std::cout << "X after 2nd stage before limit = ";
   //    std::cout << "rank = " << myRank << ", x = " << x[0] << std::endl;
   // }  
   // limit x
   
   limiter.limit(x);
 
   // if (myRank == 0)
   // {
   //    std::cout << "X after 2nd stage = ";
   //    std::cout << x[0] << std::endl;
   // }
   //    x.Print(std::cout);

   // update time step
   t += dt;

   // cout << "end step" << endl;
}

