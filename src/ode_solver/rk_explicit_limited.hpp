#ifndef RK_EXPLICIT_LIMITED_H
#define RK_EXPLICIT_LIMITED_H

#include "mfem.hpp"
#include "limiter.hpp"

using namespace std;
using namespace mfem;

/// Proc rank 
extern int myRank;

/// 
/// An explicit Runge-Kutta method corresponding to a general Butcher tableau
/// with the DG limitation
///    +--------+----------------------+
    // | c[0]   | a[0]                 |
    // | c[1]   | a[1] a[2]            |
    // | ...    |    ...               |
    // | c[s-2] | ...   a[s(s-1)/2-1]  |
    // +--------+----------------------+
    // |        | b[0] b[1] ... b[s-1] |
    // +--------+----------------------+
class ExplicitRKLimitedSolver : public ODESolver
{
private:

   /// Order
   int s;

   /// Butcher tableau coefficients
   const double *a, *b, *c;

   /// Intermediate solution on RK studies
   Vector y;

   /// Vector of rhs on different RK stages
   Vector *k;

   /// Slope limiter
   Limiter& limiter;

public:
   
   // Constructor
   ExplicitRKLimitedSolver(int _s, const double *_a, const double *_b,
                    const double *_c, Limiter& _l);

   /// Init memory
   virtual void Init(TimeDependentOperator &_f);

   /// Time step
   //  @param x previous solution
   //  @param t current time
   //  @param dt current time step
   virtual void Step(Vector &x, double &t, double &dt);

   /// Destructor
   virtual ~ExplicitRKLimitedSolver() {delete [] k;}
};

#endif // RK_EXPLICIT_LIMITED_H
