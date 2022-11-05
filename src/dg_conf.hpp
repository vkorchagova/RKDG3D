/// Configuration file

#ifndef DG_SOLVER_CONF_H
#define DG_SOLVER_CONF_H

#include <ryml.hpp>


/// --- Physics parameters --- ///

/// Specific heat ratio
const double DEFAULT_SPECIFIC_HEAT_RATIO = 1.4;

/// Individual gas constant
const double DEFAULT_GAS_CONSTANT = 287;

/// Covolume constant
const double DEFAULT_COVOLUME_CONSTANT = 0.0;


/// --- Meshing settings --- ///

/// Number of mesh uniform refinement levels in serial mode
/// (0 - do not run uniform mesh refinement)
const int DEFAULT_SERIAL_MESH_REF_LEVELS = 0;

/// Number of mesh uniform refinement levels in parallel mode (refine mesh parts on processors)
/// (0 - do not run uniform mesh refinement)
const int DEFAULT_PARALLEL_MESH_REF_LEVELS = 0;

/// Use adaptive mesh refinement
const bool DEFAULT_USE_ADAPTIVE_MESH_REFINIMENT = false;

/// Use local mesh refinement (is some mesh region)
const bool DEFAULT_USE_LOCAL_MESH_REFINEMENT = false;


/// --- Numerical parameters --- ///

/// Number of equations
const int DEFAULT_NUMBER_OF_EQUATIONS = 4;

/// Turn on restart from some previous time step
const bool DEFAULT_CHECK_FOR_RESTART = false;

/// Cut nonlinear terms when limit solution
const bool DEFAULT_LINEARIZE_SOLUTION = false;

/// Cut all slopes in case of non-physical values in cell nodes after limiting
const bool DEFAULT_REMOVE_SLOPES_AFTER_LIMITING = true;

/// Name for the empty cell group when use ony piecewise-constant solution
const ryml::csubstr DEFAULT_NAME_OF_EMPTY_GROUP_WITHOUT_SLOPES = "";

/// Small number used in near-zero thresholds and prevent division by zero
const double DEFAULT_SMALL_EPSILON = 1e-6;

/// Large number used in min/max search
const double DEFAULT_LARGE_NUMBER = 1e+9;

/// Define stencil size for preliminar memory reservation
const int DEFAULT_MAX_STENCIL_SIZE = 8;

/// Threshold value for indicator: if less we should correct slopes
const double DEFAULT_INDICATOR_CORRECTION_THRESHOLD_VALUE = 0.99999;

/// Percentage of max solution value to check diff for large values (like energy)
const double DEFAULT_BJ_DIFF_MAX_PERCENT = 1e-3;

/// Index for root component of solution to compute BJ limiter value
const int DEFAULT_BJ_ROOT_SOL_COMP_INDEX = 0;

/// Y threshold used in Michalak limiter
/// Do not used as user setting: value is chosen according to Michalak's paper.
const double DEFAULT_MICHALAK_YSTAR = 1.5;

/// Numerical constant in Shu indicator
/// Do not used as user setting: value is chosen according to Shu's paper.
const double DEFAULT_SHU_CK = 0.03;


/// --- Postprocessing options --- ///

/// Name of the folder for output results
const std::string DEFAULT_OUTPUT_FOLDER_NAME = "PV";

/// Turn on total energy computation and output
const bool DEFAULT_COMPUTE_TOTAL_ENERGY = false;

/// Write indiator values as output field
const bool DEFAULT_WRITE_INDICATORS = false;

/// Number of step intervals for visualisation
/// (1 - visu each step, 2 - visu each 2nd step, etc.)
const int DEFAULT_NUMBER_OF_OUTPUT_STEPS = 1;

/// Level of details in ParaView
/// (1 - use only cell center values, 2 - use node values ...)
const int DEFAULT_LEVEL_OF_PARAVIEW_DETAILS = 2;


#endif // DG_SOLVER_CONF_H