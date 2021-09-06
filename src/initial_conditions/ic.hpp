#ifndef IC_H
#define IC_H

#include "mfem.hpp"
#include "physics.hpp"

using namespace std;
using namespace mfem;

extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Initial conditions 
///

class IC
{

protected:

    std::function<double(const Vector& x)> initRho;
    std::function<double(const Vector& x)> initP;
    std::function<double(const Vector& x)> initU;
    std::function<double(const Vector& x)> initV;
    std::function<double(const Vector& x)> initW;

public:
    IC() {};
    virtual ~IC() {};

    virtual std::function<void(const Vector &, Vector &)> setIC() = 0;
};


class ICConstant : public IC
{
    double den;
    double velX;
    double velY;
    double velZ;
    double pres;

public:

    ICConstant(const Vector& sol) 
     : 
    IC(),
    den(sol[0]),
    velX(sol[1]),
    velY(sol[2]),
    velZ(sol.Size() == 5 ? sol[3] : 0.0),
    pres(sol[num_equation-1])
    {};
    
    ~ICConstant() {};

    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) { return den; };
        initP   = [=](const Vector& x) { return pres;  };
        initU   = [=](const Vector& x) { return velX; };
        initV   = [=](const Vector& x) { return velY; };
        initW   = [=](const Vector& x) { return velZ; };

        std::function<void(const Vector &, Vector &)> F = \
            [=](const Vector &x, Vector &y)
            {
                const int dim = x.Size();
                if (dim == 2)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        0.0, 
                        initP(x));
                }
                else if (dim == 3)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = initRho(x)*initW(x);
                    y(4) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        initW(x), 
                        initP(x));
                }
            };
    
        return F;
    }
};

class ICPlaneBreakup : public IC
{
    Vector origin;
    Vector normal;

    double den1;
    double velX1;
    double velY1;
    double velZ1;
    double pres1;

    double den2;
    double velX2;
    double velY2;
    double velZ2;
    double pres2;

public:

    ICPlaneBreakup(const Vector& sol1_, const Vector& sol2_,const Vector& origin_, const Vector& normal_) 
     : 
    IC(),
    den1(sol1_[0]),
    velX1(sol1_[1]),
    velY1(sol1_[2]),
    velZ1(sol1_.Size() == 5 ? sol1_[3] : 0.0),
    pres1(sol1_[num_equation-1]),
    den2(sol2_[0]),
    velX2(sol2_[1]),
    velY2(sol2_[2]),
    velZ2(sol2_.Size() == 5 ? sol2_[3] : 0.0),
    pres2(sol2_[num_equation-1]),
    origin(origin_),
    normal(normal_)
    {
        // origin = origin_;
        // normal = normal_;
        // origin(origin_.GetData(),origin_.Size()),
        // normal(normal_.GetData(),normal_.Size())
    };
    
    ~ICPlaneBreakup() {};

    double plane(const Vector& x)
    {
        const int dim = x.Size();
        if (dim == 2)
            return normal[0] * (x[0] - origin[0]) + \
                    normal[1] * (x[1] - origin[1]);
        else if (dim == 3)
            return normal[0] * (x[0] - origin[0]) + \
                    normal[1] * (x[1] - origin[1]) + \
                    normal[2] * (x[2] - origin[2]);
        else
        {
            cout << "Wrong dimension inside planeBreakup" << endl;
            exit(1);
        }
    }

    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) 
        { 
            return plane(x) < 0 ? den1 : den2; 
        };
        initP   = [=](const Vector& x) 
        { 
            return plane(x) < 0 ? pres1 : pres2;  
        };
        initU   = [=](const Vector& x) 
        { 
            return plane(x) < 0 ? velX1 : velX2; 
        };
        initV   = [=](const Vector& x) 
        { 
            return plane(x) < 0 ? velY1 : velY2; 
        };
        initW   = [=](const Vector& x) 
        { 
            return plane(x) < 0 ? velZ1 : velZ2; 
        };

        std::function<void(const Vector &, Vector &)> F = \
            [=](const Vector &x, Vector &y)
            {
                const int dim = x.Size();
                if (dim == 2)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        0.0, 
                        initP(x));
                }
                else if (dim == 3)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = initRho(x)*initW(x);
                    y(4) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        initW(x), 
                        initP(x));
                }
            };
    
        return F;
    }
};



class ICSphericalBreakup : public IC
{
    Vector origin;
    double rSquared;

    double den1;
    double velX1;
    double velY1;
    double velZ1;
    double pres1;

    double den2;
    double velX2;
    double velY2;
    double velZ2;
    double pres2;

public:

    ICSphericalBreakup(const Vector& sol1_, const Vector& sol2_,const Vector& origin_, const double radius_) 
     : 
    IC(),
    den1(sol1_[0]),
    velX1(sol1_[1]),
    velY1(sol1_[2]),
    velZ1(sol1_.Size() == 5 ? sol1_[3] : 0.0),
    pres1(sol1_[num_equation-1]),
    den2(sol2_[0]),
    velX2(sol2_[1]),
    velY2(sol2_[2]),
    velZ2(sol2_.Size() == 5 ? sol2_[3] : 0.0),
    pres2(sol2_[num_equation-1]),
    origin(origin_),
    rSquared(radius_*radius_)
    {
        // origin = origin_;
        // normal = normal_;
        // origin(origin_.GetData(),origin_.Size()),
        // normal(normal_.GetData(),normal_.Size())
    };
    
    ~ICSphericalBreakup() {};

    double sphere(const Vector& x)
    {
        const int dim = x.Size();
        if (dim == 2)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                    (x[1] - origin[1])*(x[1] - origin[1]) - rSquared;
        else if (dim == 3)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                   (x[1] - origin[1])*(x[1] - origin[1]) + \
                   (x[2] - origin[2])*(x[2] - origin[2]) - rSquared;
        else
        {
            cout << "Wrong dimension inside planeBreakup" << endl;
            exit(1);
        }
    }

    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? den1 : den2; 
        };
        initP   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? pres1 : pres2;  
        };
        initU   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velX1 : velX2; 
        };
        initV   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velY1 : velY2; 
        };
        initW   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velZ1 : velZ2; 
        };

        std::function<void(const Vector &, Vector &)> F = \
            [=](const Vector &x, Vector &y)
            {
                const int dim = x.Size();
                if (dim == 2)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        0.0, 
                        initP(x));
                }
                else if (dim == 3)
                {
                    y(0) = initRho(x);
                    y(1) = initRho(x)*initU(x);
                    y(2) = initRho(x)*initV(x);
                    y(3) = initRho(x)*initW(x);
                    y(4) = ComputeEnergy(
                        initRho(x), 
                        initU(x), 
                        initV(x), 
                        initW(x), 
                        initP(x));
                }
            };
    
        return F;
    }
};


// // Initial condition
// void setInitialConditions(const Vector &x, Vector &y)
// {
   
//    const int dim = x.Size();

//    double den = 0.0;
//    double velX = 0.0;
//    double velY = 0.0;
//    double velZ = 0.0;
//    double pres = 0.0;
//    double energy = 0.0;
//    double vel2 = 0.0;
//    double shrinv1 = 1.0 / (specific_heat_ratio - 1.);

//    if (problem == 1)
//    { // Circle Sod problem
      
//       double radius = 0.4;
//       const double xc = 0.0, yc = 0.0, zc = 0.0;

//       // double velX = 0.0;
//       // double velY = 0.0;
//       // double velZ = 0.0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;

//          den  = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) < radius*radius ? 1.0 : 0.125;// * (x(0) + 1.0);
//          pres = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) < radius*radius ? 1.0 : 0.1;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {  
//          vel2 = velX * velX + velY * velY + velZ * velZ;

//          den  = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) + (x(2) - zc)*(x(2) - zc) < radius*radius ? 1.0 : 0.125;
//          pres = (x(0) - xc)*(x(0) - xc) + (x(1) - yc)*(x(1) - yc) + (x(2) - zc)*(x(2) - zc) < radius*radius ? 1.0 : 0.1;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//    }
//    else if (problem == 2)
//    { // Forward Step problem
      
//       velX = 0.0;//3.0;//0.675;//1.65;
//       velY = 0.0;
//       velZ = 0.0;


//       den = 1.0;
//       pres = 1.0 / specific_heat_ratio;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {         
//          mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
//          vel2 = velX * velX + velY * velY + velZ * velZ;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//    }
//    else if (problem == 3)
//    { // Shu Osher problem
      
//       const double xc = 0.125;

      
//       velY = 0.0;
//       velZ = 0.0;


//       den  = x(0) < xc ? 3.857143 : 1 + 0.2*sin(8*2.0 * 3.14159265*x(0));// * (x(0) + 1.0);
//       pres = x(0) < xc ? 10.33333 : 1;
//       velX = x(0) < xc ? 2.629369 : 0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {         
//          mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
//          vel2 = velX * velX + velY * velY + velZ * velZ;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//    }
//    else if (problem == 4)
//    { // Double Mach problem
      
//       const double xc = 0.15;

      
//       velY = 0.0;
//       velZ = 0.0;


//       den  = x(0) < xc ? 8.0 : 1.4;
//       pres = x(0) < xc ? 116.518 : 1.0;
//       velX = x(0) < xc ? 8.25 : 0.0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {         
//          mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
//          vel2 = velX * velX + velY * velY + velZ * velZ;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//    }
//    else if (problem == 5)
//    { // Sod Covolume problem

//       covolume_constant = 0.001;
//       specific_heat_ratio = 1.3;
      
//       const double xc = 0.4;

//       // double velX = 0.0;
//       // double velY = 0.0;
//       // double velZ = 0.0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;

//          den  = x(0) < xc ? 100.0 : 1.0;// * (x(0) + 1.0);
//          pres = x(0) < xc ? 1e8 : 1e5;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {  
//          vel2 = velX * velX + velY * velY + velZ * velZ;

//          den  = x(0) < xc ? 100 : 1.0;
//          pres = x(0) < xc ? 1e8 : 1e5;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//    }
//    else if (problem == 6)
//    { // State Contact problem
      
//       const double xc = 0.5;

//       double velX = 0.0;
//       double velY = 0.0;
//       double velZ = 0.0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;

//          den  = x(0) < xc ? 1.4 : 1.0;// * (x(0) + 1.0);
//          pres = 1.0;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {  
//          vel2 = velX * velX + velY * velY + velZ * velZ;

//          den  = x(0) < xc ? 1.4 : 1.0;// * (x(0) + 1.0);
//          pres = 1.0;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//    }
//    else if (problem == 7)
//    { // Strip Sod problem
      
//       double stripLen = 1.0;
//       const double xc = 0.0;

//       // double velX = 0.0;
//       // double velY = 0.0;
//       // double velZ = 0.0;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;

//          den  = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.125;// * (x(0) + 1.0);
//          pres = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.1;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {  
//          vel2 = velX * velX + velY * velY + velZ * velZ;

//          den  = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.125;
//          pres = fabs(x(0) - xc) < 0.5 * stripLen ? 1.0 : 0.1;
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//    }
//    else if (problem == 8)
//    { // AstroTest

//       specific_heat_ratio = 5.0/3.0;

//       double u = 1.02; // radial velocity
//       double kTilde = 5.5e6;
//       double M0 = 1e30;
//       double G0 = 6.67e-11;
//       double Ror = 696e6;

//       double k = kTilde * pow(M0, specific_heat_ratio - 2.0) / G0 / pow(Ror, 3.0*specific_heat_ratio - 4.0);
//       double frackgamma = (specific_heat_ratio - 1.0) / k / specific_heat_ratio;

//         // source = [=](const numvector<double, dimPh> sol, const Point& r)
//         // {
//         //     return numvector<double, 5> \
//         //     { \
//         //         0.0, \
//         //         - x(0) / pow(x.Norml2(),3), \
//         //         - x(1) / pow(x.Norml2(),3), \
//         //         0.0, \
//         //         0.0 \
//         //     }; 
//         // };

//       if (dim == 2)
//       {
//          velX = - u * x(1) / pow( x.Norml2() * x.Norml2(), 0.75);
//          velY = u * x(0) / pow( x.Norml2() * x.Norml2(), 0.75);

//          vel2 = velX * velX + velY * velY;

//          den  = pow( frackgamma * (u*u - 1.0) * (x.Norml2() - 1.0) / x.Norml2(), 1.0 / (specific_heat_ratio - 1.0));
//          pres = k * pow( frackgamma * (u*u - 1.0) * (x.Norml2() - 1.0) / x.Norml2(), specific_heat_ratio / (specific_heat_ratio - 1.0));
//          energy = shrinv1 * pres / den * (1 - den * covolume_constant) + 0.5 * vel2;
//       }
//       else
//       {
//          cout << "wrong dim for AstroTest" << endl;
//       }
//    }
//    else if (problem == 9)
//    { // Forward Step problem
      
//       velX = 577.2273;//0.675;//1.65;
//       velY = 0.0;
//       velZ = 0.0;

//       den = 1.176413;
//       pres = 101325;

//       if (dim == 2)
//       {
//          vel2 = velX * velX + velY * velY;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//       else if (dim == 3)
//       {         
//          mfem_error("Cannot initialize 3D problem for ForwardStep conditions");
//          vel2 = velX * velX + velY * velY + velZ * velZ;
//          energy = shrinv1 * pres / den + 0.5 * vel2;
//       }
//    }
//    else
//    {
//       mfem_error("Cannot recognize problem."
//                  "Options are: 1 - Circle Sod problem, 2 - ForwardStep, 3 - Shu-Osher problem, 4 - DoubleMach problem, 5 - Covolume Sod problem, 6 - Contact State problem, 7 - StripSod, 8 - AstroTest, 9 - Wing");
//    }   

//    // set conservative variables
//    if (dim == 2)
//    {
//       y(0) = den;
//       y(1) = den * velX;
//       y(2) = den * velY;
//       y(3) = den * energy;
//    }
//    else if (dim == 3)
//    {
//       y(0) = den;
//       y(1) = den * velX;
//       y(2) = den * velY;
//       y(3) = den * velZ;
//       y(4) = den * energy;
//    }
// }

#endif // IC_H