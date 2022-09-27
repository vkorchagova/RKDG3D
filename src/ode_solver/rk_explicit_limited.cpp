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
   for (int i = 0; i < y.Size(); ++i) 
   {
      y(i) = 0.0;
   }
   for (int i = 0; i < s; i++)
   {
      k[i].SetSize(n, mem_type);
   }
}


void ExplicitRKLimitedSolver::Step(Vector &x, double &t, double &dt)
{
   for (int i = 0; i < s; i++)
   {
      k[i].SetSize(x.Size(), mem_type);
   }

   y.SetSize(x.Size(), mem_type);
   for (int i = 0; i < y.Size(); ++i) 
   {
      y(i) = 0.0;
   }

   f->SetTime(t);
   f->Mult(x, k[0]);
   
   for (int l = 0, i = 1; i < s; i++)
   {
      add(x, a[l++]*dt, k[0], y);
      
      for (int j = 1; j  < i; j++)
      {
         y.Add(a[l++]*dt, k[j]);
      }
      
      // limit y and compute rhs with good y
      limiter.update(y);
      
      f->SetTime(t + c[i-1]*dt);
      f->Mult(y, k[i]);
   }

   for (int i = 0; i < s; i++)
   {
      x.Add(b[i]*dt, k[i]);
   }
   
   limiter.update(x);
}

