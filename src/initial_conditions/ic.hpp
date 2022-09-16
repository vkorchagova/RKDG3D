#ifndef IC_H
#define IC_H

#include "mfem.hpp"
#include "physics.hpp"

using namespace std;
using namespace mfem;

/// Number of equations
extern int num_equation;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for different types of initial conditions 
///
class IC
{

protected:

    /// Density function
    std::function<double(const Vector& x)> initRho;

    /// Pressure function
    std::function<double(const Vector& x)> initP;

    /// velocity X component function
    std::function<double(const Vector& x)> initU;

    /// velocity Y component function
    std::function<double(const Vector& x)> initV;

    /// velocity Z component function
    std::function<double(const Vector& x)> initW;

public:

    /// Constructor
    IC() {};

    /// Destructor
    virtual ~IC() {};

    /// Define initial conditions
    virtual std::function<void(const Vector &, Vector &)> setIC() = 0;
};

/// 
/// Constant initial values in the whole flow domain
///
class ICConstant : public IC
{
    /// Initial value of density
    double den;

    /// Initial value of velocity X component
    double velX;

    /// Initial value of velocity Y component
    double velY;

    /// Initial value of velocity Z component
    double velZ;

    /// Initial value of pressure
    double pres;

public:

    /// Constructor
    ICConstant(const Vector& sol) 
     : 
    IC(),
    den(sol[0]),
    velX(sol[1]),
    velY(sol[2]),
    velZ(sol.Size() == 5 ? sol[3] : 0.0),
    pres(sol[num_equation-1])
    {};
    
    /// Destructor
    ~ICConstant() {};

    /// Define initial conditions
    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) { return den; };
        initP   = [=](const Vector& x) { return pres; };
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

/// 
/// Two constant initial states in the flow domain separated by plane
///
class ICPlaneBreakup : public IC
{
    /// Plane origin
    Vector origin;

    /// Plane normal
    Vector normal;

    /// Initial value of density behind plane
    double den1;

    /// Initial value of velocity X component behind plane
    double velX1;

    /// Initial value of velocity Y component behind plane
    double velY1;

    /// Initial value of velocity Z component behind plane
    double velZ1;

    /// Initial value of pressure behind plane
    double pres1;

    /// Initial value of density forward of plane
    double den2;

    /// Initial value of velocity X component forward of plane
    double velX2;

    /// Initial value of velocity Y component forward of plane
    double velY2;

    /// Initial value of velocity Z component forward of plane
    double velZ2;

    /// Initial value of pressure forward of plane
    double pres2;

public:

    /// Constructor
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
    {};
    
    /// Destructor
    ~ICPlaneBreakup() {};

    /// Define a plane equation
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
            cout << "Wrong dimension inside ICPlaneBreakup" << endl;
            exit(1);
        }
    }

    /// Define initial conditions
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


/// 
/// Two constant initial states in the flow domain separated by sphere
///
class ICSphericalBreakup : public IC
{
    /// Center of sphere
    Vector origin;

    /// Radius of sphere (squared)
    double rSquared;

    /// Radius of smoothing part (linear transition between two states) (squared)
    double rSquaredDelta;

    /// Initial value of density inside sphere
    double den1;

    /// Initial value of velocity X component inside sphere
    double velX1;

    /// Initial value of velocity Y component inside sphere
    double velY1;

    /// Initial value of velocity Z component inside sphere
    double velZ1;

    /// Initial value of pressure inside sphere
    double pres1;

    /// Initial value of density outside sphere
    double den2;

    /// Initial value of velocity X component outside sphere
    double velX2;

    /// Initial value of velocity Y component outside sphere
    double velY2;

    /// Initial value of velocity Z component outside sphere
    double velZ2;

    /// Initial value of pressure outside sphere
    double pres2;

public:

    /// Constructor
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
    rSquared(radius_*radius_),
    rSquaredDelta(0.2*0.2*rSquared)
    {};
    
    /// Destructor
    ~ICSphericalBreakup() {};

    /// Define a sphere equation
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
            cout << "Wrong dimension inside ICSphericalBreakup" << endl;
            exit(1);
        }
    }

    /// Get a linear transition between two states
    double linear(const Vector& x, double r1, double r2, double val1, double val2)
    {
        return (x.DistanceTo(origin) - r1) * (val2 - val1) / (r2 - r1) + val1;
    }

    /// Define initial conditions
    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? den1 : (sphere(x) < rSquaredDelta ? linear(x, sqrt(rSquared), sqrt(rSquared) + sqrt(rSquaredDelta), den1, den2) : den2); 
        };
        initP   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? pres1 : (sphere(x) < rSquaredDelta ? linear(x, sqrt(rSquared), sqrt(rSquared) + sqrt(rSquaredDelta), pres1, pres2) : pres2);  
        };
        initU   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velX1 : (sphere(x) < rSquaredDelta ? linear(x, sqrt(rSquared), sqrt(rSquared) + sqrt(rSquaredDelta), velX1, velX2) : velX2); 
        };
        initV   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velY1 : (sphere(x) < rSquaredDelta ? linear(x, sqrt(rSquared), sqrt(rSquared) + sqrt(rSquaredDelta), velY1, velY2) : velY2); 
        };
        initW   = [=](const Vector& x) 
        { 
            return sphere(x) < 0 ? velZ1 : (sphere(x) < rSquaredDelta ? linear(x, sqrt(rSquared), sqrt(rSquared) + sqrt(rSquaredDelta), velZ1, velZ2) : velZ2); 
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

/// 
/// Constant initial values in the whole flow domain for velocity and pressure
/// Exponential pulse for density
///
class ICDensityPulse : public IC
{
    /// Origin of density pulse
    Vector origin;

    /// Amplitude of density pulse
    double epsilon;

    /// Initial value of velocity X component
    double velX;

    /// Initial value of velocity Y component
    double velY;

    /// Initial value of velocity Z component
    double velZ;

    /// Initial value of pressure
    double pres;

    /// Compute 1-r^2 for defined point
    double oneMinusR2(const Vector& x)
    {
        const int dim = x.Size();
        if (dim == 2)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                    (x[1] - origin[1])*(x[1] - origin[1]) - 1.0;
        else if (dim == 3)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                   (x[1] - origin[1])*(x[1] - origin[1]) + \
                   (x[2] - origin[2])*(x[2] - origin[2]) - 1.0;
        else
        {
            cout << "Wrong dimension inside ICDensityPulse" << endl;
            exit(1);
        }
    }

public:

    ICDensityPulse(const Vector& sol_, const Vector& origin_, const double epsilon_) 
     : 
    IC(),
    velX(sol_[1]),
    velY(sol_[2]),
    velZ(sol_.Size() == 5 ? sol_[3] : 0.0),
    pres(sol_[num_equation-1]),
    origin(origin_),
    epsilon(epsilon_)
    {};
    
    ~ICDensityPulse() {};

    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) 
        { 
            return pow((1.0 - (specific_heat_ratio - 1.0)*epsilon*epsilon*exp(oneMinusR2(x))/(8.0*specific_heat_ratio*3.14169265*3.14169265)), 1.0 / (specific_heat_ratio - 1.0)); 
        };
        initP   = [=](const Vector& x) 
        { 
            return pres;  
        };
        initU   = [=](const Vector& x) 
        { 
            return velX; 
        };
        initV   = [=](const Vector& x) 
        { 
            return velY; 
        };
        initW   = [=](const Vector& x) 
        { 
            return velZ; 
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



/// 
/// Two constant initial states in the flow domain separated by two hexahedrons
///
class ICTwoHexahedronsBreakup : public IC
{
    /// min point, first hex
    Vector min1;

    /// max point, first hex
    Vector max1;

    /// min point, second hex
    Vector min2;

    /// max point, second hex
    Vector max2;

    /// Initial value of density inside hexahedrons
    double den1;

    /// Initial value of velocity X component inside hexahedrons
    double velX1;

    /// Initial value of velocity Y component inside hexahedrons
    double velY1;

    /// Initial value of velocity Z component inside hexahedrons
    double velZ1;

    /// Initial value of pressure inside hexahedrons
    double pres1;

    /// Initial value of density outside hexahedrons
    double den2;

    /// Initial value of velocity X component outside hexahedrons
    double velX2;

    /// Initial value of velocity Y component outside hexahedrons
    double velY2;

    /// Initial value of velocity Z component outside hexahedrons
    double velZ2;

    /// Initial value of pressure outside hexahedrons
    double pres2;

public:

    /// Constructor
    ICTwoHexahedronsBreakup(const Vector& sol1_, const Vector& sol2_, const Vector& min1_, const Vector& max1_, const Vector& min2_, const Vector& max2_) 
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
    min1(min1_),
    max1(max1_),
    min2(min2_),
    max2(max2_)
    {};
    
    /// Destructor
    ~ICTwoHexahedronsBreakup() {};

    /// Define region inside hexadrons
    bool inHexadrons(const Vector& x)
    {
        const int dim = x.Size();
        if (dim == 2)
            return ((x[0] < max1[0] && x[0] > min1[0]) &&
                    (x[1] < max1[1] && x[1] > min1[1]) ) ||
                   ((x[0] < max2[0] && x[0] > min2[0]) &&
                    (x[1] < max2[1] && x[1] > min2[1]) );
        else if (dim == 3)
            return ((x[0] < max1[0] && x[0] > min1[0]) &&
                    (x[1] < max1[1] && x[1] > min1[1]) &&
                    (x[2] < max1[2] && x[2] > min1[2]) ) ||
                   ((x[0] < max2[0] && x[0] > min2[0]) &&
                    (x[1] < max2[1] && x[1] > min2[1]) &&
                    (x[2] < max2[2] && x[2] > min2[2]) );
        else
        {
            cout << "Wrong dimension inside ICTwoHexahedronsBreakup" << endl;
            exit(1);
        }
    }

    /// Define initial conditions
    virtual std::function<void(const Vector &, Vector &)> setIC() override
    {
        initRho = [=](const Vector& x) 
        { 
            return inHexadrons(x) ? den1 : den2; 
        };
        initP   = [=](const Vector& x) 
        { 
            return inHexadrons(x) ? pres1 : pres2;  
        };
        initU   = [=](const Vector& x) 
        { 
            return inHexadrons(x) ? velX1 : velX2; 
        };
        initV   = [=](const Vector& x) 
        { 
            return inHexadrons(x) ? velY1 : velY2; 
        };
        initW   = [=](const Vector& x) 
        { 
            return inHexadrons(x) ? velZ1 : velZ2; 
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




#endif // IC_H