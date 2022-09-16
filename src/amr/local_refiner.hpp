#ifndef LOCAL_REFINER_H
#define LOCAL_REFINER_H
 
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/// Proc rank 
extern int myRank;

/// 
/// Abstract class for different types of local refiner 
///
class LocalRefiner
{

protected:

    /// Pointer to mesh
    Mesh* mesh;

    /// Refinement depth
    int level_of_refinement;

    /// Refine inside or outside
    bool ref_inside;

    /// Marked cells to be refined
    Array<int> marked_cells;

    /// Element vertices
    Array<int> el_vertices;

    /// Space dimension
    int dim;

public:

    /// Constructor
    LocalRefiner
    (
        Mesh* _mesh, 
        int _lor, 
        bool _in
    ) : 
        mesh(_mesh), 
        level_of_refinement(_lor), 
        ref_inside(_in),
        dim(mesh->SpaceDimension())
    {};

    /// Destructor
    virtual ~LocalRefiner() { marked_cells.DeleteAll(); };

    /// Refinement domain
    virtual bool inside_domain(const double* x) = 0;

    /// Mark mesh cells to be refined
    void MarkCells()
    {
        for (int iCell = 0; iCell < mesh->GetNE(); iCell++)
        {   
            mesh->GetElementVertices(iCell, el_vertices);

            for (int iVert : el_vertices)
                if ( !(inside_domain(mesh->GetVertex(iVert)) ^ ref_inside) )
                {
                    marked_cells.Append(iCell);
                    break;
                }
        }
        // marked_cells.Print (cout);
    }

    /// Refinement process
    void Refine()
    {
        for (int l = 0; l < level_of_refinement; l++)
        {
            marked_cells.DeleteAll();
            MarkCells();
            mesh->GeneralRefinement(marked_cells, 1, 1); // 1 = non-conforming mesh, 1 = only one hanging node
        }
    }
};

/// 
/// Refine mesh cells inside (or outside) a sphere
///
class SphericalLocalRefiner : public LocalRefiner
{
    /// Center of sphere
    Vector origin;

    /// Radius of sphere (squared)
    double rSquared;

public:

    /// Constructor
    SphericalLocalRefiner
    (
        const Vector& origin_, 
        const double radius_,
        Mesh* _mesh, 
        int _lor, 
        bool _in
    ) : 
        LocalRefiner(_mesh, _lor, _in),
        origin(origin_),
        rSquared(radius_*radius_)
    {};
    
    /// Destructor
    virtual ~SphericalLocalRefiner() {};

    /// Define a sphere equation
    virtual bool inside_domain(const double* x) override
    {
        if (dim == 2)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                    (x[1] - origin[1])*(x[1] - origin[1]) - rSquared - 1e-6 < 0.0;
        else if (dim == 3)
            return (x[0] - origin[0])*(x[0] - origin[0]) + \
                   (x[1] - origin[1])*(x[1] - origin[1]) + \
                   (x[2] - origin[2])*(x[2] - origin[2]) - rSquared - 1e-6 < 0.0;
        else
        {
            cout << "Wrong dimension inside SphericalLocalRefiner" << endl;
            exit(1);
        }
    }
};


/// 
/// Refine mesh cells inside (or outside) a hexahedron
///
class HexahedronLocalRefiner : public LocalRefiner
{
    /// min point
    Vector min1;

    /// max point
    Vector max1;

public:

    /// Constructor
    HexahedronLocalRefiner
    (
        const Vector& min1_, 
        const Vector& max1_,
        Mesh* _mesh, 
        int _lor, 
        bool _in
    ) : 
        LocalRefiner(_mesh, _lor, _in),
        min1(min1_),
        max1(max1_)
    {};
    
    /// Destructor
    virtual ~HexahedronLocalRefiner() {};

    /// Define region inside hexadron
    virtual bool inside_domain(const double* x) override
    {
        if (dim == 2)
            return ((x[0] < max1[0] && x[0] > min1[0]) &&
                    (x[1] < max1[1] && x[1] > min1[1]) );
        else if (dim == 3)
            return ((x[0] < max1[0] && x[0] > min1[0]) &&
                    (x[1] < max1[1] && x[1] > min1[1]) &&
                    (x[2] < max1[2] && x[2] > min1[2]) );
        else
        {
            cout << "Wrong dimension inside HexahedronLocalRefiner" << endl;
            exit(1);
        }
    }
};


/// 
/// Refine mesh cells inside (or outside) a cylinder
///
class CylindricalLocalRefiner : public LocalRefiner
{
    /// min point
    Vector p1;

    /// max point
    Vector p2;

    /// radius*radius
    double rSquared;

public:

    /// Constructor
    CylindricalLocalRefiner
    (
        const Vector& p1_, 
        const Vector& p2_,
        double radius_,
        Mesh* _mesh, 
        int _lor, 
        bool _in
    ) : 
        LocalRefiner(_mesh, _lor, _in),
        p1(p1_),
        p2(p2_),
        rSquared(radius_*radius_)
    {};
    
    /// Destructor
    virtual ~CylindricalLocalRefiner() {};

    /// Define region inside hexadron
    virtual bool inside_domain(const double* x) override
    {
        Vector p(3);
        subtract(p2,p1,p);
        Vector xv(3); 
        xv[0] = x[0];
        xv[1] = x[1];
        xv[2] = x[2];
        Vector a(3);
        subtract(xv,p1,a);
        double magAP = a*p;
        double magP1P2Sq = p2.DistanceSquaredTo(p1);
        if (magAP > 0 && magAP < magP1P2Sq)
        {
            double d = xv.DistanceSquaredTo(p1) - magAP*magAP / magP1P2Sq;
            return d < rSquared;
        }
        return 0;
    }
};

#endif // LOCAL_REFINER_H