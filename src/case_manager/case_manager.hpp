#ifndef CASE_MANAGER_H
#define CASE_MANAGER_H

#include "mfem.hpp"
#include "ic.hpp"
#include "rk_explicit_limited.hpp"

#include "rs_rusanov.hpp"
#include "rs_hll.hpp"
#include "rs_hllc.hpp"
#include "rs_llf.hpp"

#include "boundary_integrator_wall.hpp"
#include "boundary_integrator_char_outlet.hpp"
#include "boundary_integrator_open.hpp"
#include "boundary_integrator_open_fixed_pressure.hpp"
#include "boundary_integrator_open_total_pressure.hpp"
#include "boundary_integrator_subsonic_inlet_total_pressure.hpp"
#include "boundary_integrator_subsonic_inlet_fixed_pressure.hpp"
#include "boundary_integrator_supersonic_inlet.hpp"

#include "limiter_findiff.hpp"
#include "limiter_multiplier.hpp"
#include "indicator_nowhere.hpp"
#include "indicator_everywhere.hpp"
#include "indicator_bj.hpp"
#include "indicator_shu.hpp"
#include "averager.hpp"

#include <filesystem>
#include <queue>
#include <map>
#include <ryml_std.hpp>
#include <ryml.hpp>

using namespace mfem;

extern int problem;
extern int num_equation;
extern int myRank;

extern double specific_heat_ratio;
extern double covolume_constant;
extern double gas_constant;

/// 
/// Case manager class 
/// Reads case YAML file, initialize parameters, run computations
///

class CaseManager
{
   /// Pointer to full content of YAML file
   std::string* contents;

   /// True if needed restart from non-zero time
   bool restart;

   /// Time step for restart
   int restartCycle;

   /// Window for saved frames 
   int nSavedFrames;

   /// VisIt data collection for restart
   VisItDataCollection& restart_data_c;

   /// Queue for restart time step names
   std::priority_queue<std::string, std::vector<std::string>, std::greater<std::string>> restart_queue;

   /// Restart names for window moving
   std::string restart_current_cycle_name;
   std::string restart_deleted_cycle_name;
   // stream for formation of file names
   std::ostringstream restart_name_stream;

   /// Serial refinement mesh levels
   int serRefLevels;

   /// Parallel refinement mesh levels (on each processor locally)
   int parRefLevels;

   /// True if adaptively refinement mesh
   bool adaptive_mesh;

   /// Spatial polynomial order
   int spatialOrder;

   // for adaptive mesh
   double total_error_fraction;
   double max_elem_error;
   double hysteresis;
   int nc_limit ;
   bool prefer_conforming_refinement;

   // Type of initial condition
   std::string icType;

   // Type of Riemann solver
   std::string rsolverType;

   // Pointer to initial conditions interface
   IC* ICInterface;

   // Mapping of names for boundary groups and physical groups of elements
   std::map<std::string,int> map_phys_names_tag_1D;
   std::map<std::string,int> map_phys_names_tag_2D;
   std::map<std::string,int> map_phys_names_tag_3D;

   // Array of boundary markers for BdrIntegrators
   std::vector<Array<int>> bdr_markers;

   Vector origin; //(num_equation-2);
   Vector normal; //(num_equation-2)

//   double gamma;
   
   // Runge -- Kutta Butcher coefficients
   double* a;
   double* b;
   double* c;

   /// Postprocessing features
   bool checkTotalEnergy;
   bool writeIndicators;

   /// Limiter features
   bool linearize;
   bool haveLastHope;

   // Set initial conditions
   static void setIC(const Vector&x, Vector& y);

   // GMSH sometimes makes strange big numbers for attributes (for ex. 167615451)
   // MFEM creates arrays for boundary markers, their size is equal to maximal attribute value.
   // To avoid memory overload, attributes should be changed for small numbers.
   void minimizeAttributes(ParMesh*& pmesh);

   // Additional reading for physical groups for possibility to write boundaru names in setttings instead of boundary tags directly
   void readPhysicalNames(ParMesh*& mesh);

public:

   // Case directory 
   std::string caseDir;

   // YAML tree of settings
   ryml::Tree* settings;

   // Parse YAML settings file
   void parse(std::string& caseFileName);

public:

   /// Constructor
   CaseManager(std::string& caseFileName, VisItDataCollection& rdc); 

   /// Destructor
   ~CaseManager() 
   {
      delete[] c;
      delete[] b;
      delete[] a;
      delete ICInterface;
      delete settings;
      delete contents;
      origin.SetSize(1);
      origin[0] = 0;
      normal.SetSize(1);
      normal[0] = 0;
   };

   /// Load mesh from file or from previous saved data
   void loadMesh(ParMesh*& pmesh);

   /// Load adaptive mesh settings
   void loadAdaptiveMeshSettings(ThresholdRefiner& refiner,ThresholdDerefiner& derefiner);

   /// Load initial solution from function or from previous saved data
   void loadInitialSolution(ParFiniteElementSpace& vfes, const Array<int>& offsets, BlockVector& u_block, ParGridFunction& sol);

   /// Load boundary conditions
   void addBoundaryIntegrators(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim);

   /// Return spatial DG order
   int getSpatialOrder() {return spatialOrder;};

   /// Check number of equations
   void checkNumEqn(int dim);

   /// Load Riemann solver object of appropriate type
   void loadRiemannSolver(RiemannSolver*& rsolver);

   // Load indicator and limiter
   void loadLimiter(Averager& avgr, Indicator*& ind, Limiter*& l, const Array<int>& offsets, const int dim, BlockVector& indicatorData, ParFiniteElementSpace& vfes);

   // Load RK time solver
   void loadTimeSolver(ODESolver*& ode_solver, Limiter* l);

   // Load time control parameters
   void loadTimeSettings(double& t_final, double& dt, double& cfl);

   // Initialize restart queue
   void initializeRestartQueue();

   // Set start time for time cycle
   double setStartTime();
   int setStartTimeCycle();

   void cleanPreviousRestartFrames(int& ti);

   bool is_restart() {return restart;};

   bool is_adaptive() {return adaptive_mesh;};

   bool check_total_energy() {return checkTotalEnergy;};

   bool write_indicators() {return writeIndicators; }

   void getVisSteps(int& vis_steps);
};



#endif // CASE_MANAGER_H

