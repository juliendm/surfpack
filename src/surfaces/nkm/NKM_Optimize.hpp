#ifndef __OPTIMIZE_HPP__ 
#define __OPTIMIZE_HPP__ 

#include "NKM_SurfData.hpp"

namespace nkm {

class SurfPackModel;


/***********************************************************/
/**** definition of optimizer specific data starts here ****/
/***********************************************************/

/// settings specific to guess and check strategy
struct OptGuessData {
  /// number of recomended random guesses for randomly guessing parameter sets
  unsigned numGuesses;
};


///contains variables specific to CONMIN optimizer
struct OptProbConminData
{
/*  The following was copied from the conmin user's manual found at 
    http://www.eng.buffalo.edu/Research/MODEL/mdo.test.orig/CONMIN/manual.html#Param_Main
    
IPRINT   Print control.  All printing is done on unit number 6.

         0:  Print nothing.

         1:  Print initial and final function information.

         2:  1st debug level.  Print all of above plus control
             parameters.  Print function value and X-vector at each
             iteration.

         3:  2nd. debug level.  Print all of above plus all constraint
             values, numbers of active or violated constraints, direction
             vectors, move parameters and miscellaneous information.  The
             constraint parameter, BETA, printed under this option
             approaches zero as the optimum objective is achieved.

         4:  Complete debug.  Print all of above plus gradients of
             objective function, active or violated constraint functions
             and miscellaneous information.

ITMAX    Default value = 10.  Maximum number of iterations in the
         minimization process.  If NFDG.EQ.0 each iteration requires one
         set of gradient computations (INFO = 3 or 4) and approximately
         three function evaluations (INFO = 1 or 2).  If NFDG.GT.0
         each iteration requires approximately NDV + 3 function
         evaluations (INFO = 1 or 2).

NSIDE    Side constraint parameter.  NSIDE = 0 signifies that the
         variables X(I) do not have lower or upper bounds.  NSIDE.GT.0
         signifies that all variables X(I) have lower and upper bounds
         defined by VLB(I) and VUB(I) respectively.  If one or more
         variables are not bounded while others are, the values of the
         lower and upper bounds on the unbounded variables must be taken
         as very large negative and positive values respectively
         (i.e., VLB(I) = -1.0E+10, VUB(I) = 1.0E+10).

ICNDIR   Default value = NDV + 1.  Conjugate direction restart parameter.
         If the function is currently unconstrained, (all G(J).LT.CT or
         NCON = NSIDE = 0), Fletcher-Reeves conjugate direction method will
         be restarted with a steepest descent direction every ICNDIR
         iterations.  If ICNDIR = 1 only steepest descent will be used.

NFDG     Gradient calculation control parameter. 

         NFDG = 0:  all gradient information is calculated by finite difference
                    within CONMIN.

         NFDG = 1:  all gradient information is supplied by the user.

         NFDG = 2:  the gradient of OBJ is supplied by the user and the
                    gradients of constraints are calculated by finite
                    difference within CONMIN.

	 During execution, the arrays g, a, and isc must be calculated by the 
	 user, depending on the value of NFDG. If NFDG = 0, only array g is 
	 calculated. If NFDG = 2, only arrays g and df are calculated. The 
	 remaining arrays are used internally by CONMIN.

FDCH     Default value = 0.01.  Not used if NFDG = 0.  Relative change in
         decision variable X(I) in calculating finite difference
         gradients.  For example, FDCH = 0.01 corresponds to a finite
         difference step of one percent of the value of the decision
         variable.

FDCHM    Default value = 0.01.  Not used if NFDG = 0.  Minimum absolute
         step in finite difference gradient calculations.  FDCHM applies
         to the unscaled variable values.

CT       Default value = -0.1.  Not used if NCON = NSIDE = 0.
         Constraint thickness parameter.  If CT.LE.G(J).LE.ABS(CT),
         G(J) is defined as active.  If G(J).GT.ABS(CT), G(J) is said to
         be violated.  If G(J).LT.CT, G(J) is not active.  CT is
         sequentially reduced in magnitude during the optimization
         process.  If ABS(CT) is very small, one or more constraints
         may be active on one iteration and inactive on the next,
         only to become active again on a subsequent iteration.
         This is often referred to as "zigzagging" between constraints.
         A wide initial value of the constraint thickness is desirable
         for highly nonlinear problems so that when a constraint
         becomes active it tends to remain active, thus reducing the
         zigzagging problem.  The default value is usually adequate.

CTMIN    Default value = 0.004.  Not used if NCON = NSIDE = 0.  Minimum
         absolute value of CT considered in the optimization process.
         CTMIN may be considered as "numerical zero" since it may not be
         meaningful to compare numbers smaller than CTMIN.  The value of
         CTMIN is chosen to indicate that satisfaction of a constraint
         within this tolerance is acceptable.  The default value is usually
         adequate.

CTL      Default value = -0.01.  Not used if NCON = NSIDE = 0.
         Constraint thickness parameter for linear and side constraints.
         CTL is smaller in magnitude than CT because the zigzagging
         problem is avoided with linear and side constraints.  The default
         value is usually adequate.

CTLMIN   Default value = 0.001.  Not used if NCON = NSIDE = 0.  Minimum
         absolute value of CTL considered in the optimization process.
         The default value is usually adequate.

DELFUN   Default value = 0.001.  Minimum relative change in the objective
         function to indicate convergence.  If in ITRM consecutive
         iterations, ABS(1.0-OBJ(J-1)/OBJ(J)).LT.DELFUN and the current
         design is feasible (all G(J).LE.ABS(CT)), the minimization
         process is terminated.  If the current design is infeasible
         (some G(J).GT.ABS(CT)), five iterations are required to
         terminate and this situation indicates that a feasible design
         may not exist.

DABFUN   Default value = 0.001 times the initial function value.  Same
         as DELFUN except comparison is on absolute change in the
         objective function, ABS(OBJ(J)-OBJ(J-1)), instead of relative
         change.

ITRM     Default value = 3.  Number of consecutive iterations to indicate
         convergence by relative or absolute changes, DELFUN or DABFUN.


We may want to include additional CONMIN variables here
*/

  ///Finite difference flag.
  int nfdg; //=0;

  ///Flag to control amount of output data.
  int iprint;

  ///Flag to specify the maximum number of iterations.
  int itmax;

  ///Relative finite difference step size.
  double fdch; //= 1.0e-2;

  ///Absolute finite difference step size.
  double fdchm; //= 1.0e-2;

  ///Constraint thickness parameter: The value of CT decreases in magnitude during optimization
  double ct; //= -0.1;

  ///Minimum absolute value of CT used during optimization.
  double ctmin; //= 0.004;

  ///Constraint thickness parameter for linear and side constraints.
  double ctl; //= -0.01;

  ///Minimum value of CTL used during optimization.
  double ctlmin; //= 0.001;

  ///Relative convergence criterion threshold. Threshold for the minimum relative change in the objective function.
  double delfun; //= .001;

  ///Absolute convergence criterion threshold. Threshold for the minimum relative change in the objective function.
  double dabfun; //= 1.0e-10;

  ///Internal CONMIN variable: side constraints parameter
  int nside; //=1;

  ///Internal CONMIN variable: diminishing return criterion iteration number.
  int  itrm; //= 3;

  /// Internal CONMIN variable: conjugate direction restart parameter.
  int  icndir;  //=NDV+1;
};


struct OptProbDirectData {

  bool constraintsPresent;

  bool verboseOutput;

  double minBoxSize; // -1 will yield default of 1.0e-4;

  double volBoxSize; // -1 will yield default of 1.0e-6;

  double solutionTarget; // = -DBL_MAX;

  double convergenceTol; // = 1.0e-4;

  int maxFunctionEvals; // = 10000;

  int maxIterations; // = 1000;

};


/**
   definition of OptimizationProblem base class
*/
class OptimizationProblem
{

public:

  OptimizationProblem(SurfPackModel& model, int num_vars, 
		      int num_constraints = 0)
    :theModel(model), numDesignVar(num_vars), numConFunc(num_constraints)
  { 
    //printf("calling the OptProb constructor num_vars=%d numDesignVar=%d\n",
    //num_vars,numDesignVar); fflush(stdout);
    lowerBounds.newSize(numDesignVar,1);
    upperBounds.newSize(numDesignVar,1);
    initialIterates.newSize(numDesignVar,1);
    bestVars.newSize(numDesignVar,1);
  }
  
  // init functions

  /// set a single lower bound
  void lower_bound(int i, double lb);

  /// set a single upper bound
  void upper_bound(int i, double ub);
    
  /// set a single initial iterate
  void initial_iterate(int i, double x0); 

  void add_initial_iterates(MtxDbl& init_iterates_to_add);


  // run functions
  void conmin_optimize();

  // minimize unconstrained (or penalized) objective with DiRECT
  void direct_optimize();

  void best_guess_optimize(int num_guesses);

  void multistart_conmin_optimize(int num_guesses);


  // post functions

  /// retrieve best point from optimization
  const MtxDbl& best_point() const;


  /// controls for CONMIN (consider private)
  OptProbConminData conminData;

  // TODO: controls for guessing
  OptGuessData guessData;

  // controls for DiRECT
  OptProbDirectData directData;

  // return a random guess in [lowerBounds, upperBounds]
  // TODO: allow models to override?  Not sure why one would need to...
  void getRandGuess(MtxDbl& guess) const;

  /// get either an specified initial iterate or randomly generated
  /// one from list
  void retrieve_initial_iterate(int index, MtxDbl& iterate);

private:

  // helper functions
  

  // underlying optimizer implementations

  void optimize_with_conmin(MtxDbl& initial_iterate, double& final_val);

  void optimize_with_direct(double& final_val);

  /// 'fep' in Griffin-modified NCSUDirect: computes the value of the
  /// objective function (potentially at multiple points, passed by function
  /// pointer to NCSUDirect).  Include unscaling from DIRECT.
  static int 
  direct_objective_eval(int *n, double c[], double l[], double u[],
			int point[], int *maxI, int *start, int *maxfunc,
			double fvec[], int iidata[], int *iisize,
			double ddata[], int *idsize, char cdata[],
			int *icsize);
  
  static OptimizationProblem* optimizationProblemInstance;

  // data

  /// the model 
  /// (TODO: move to Command paradigm so we need not optimize on model)
  // TODO: generalize to Model&
  SurfPackModel& theModel;

  /// number of design variables
  int numDesignVar;

  /// number of constraint FUNCTIONS (does not include side
  /// constraints) NCON may be zero;
  int numConFunc; 

  /// lower bound constraints
  MtxDbl lowerBounds;

  /// upper bound constraints
  MtxDbl upperBounds;

  /// matrix of (possibly multiple) initial guesses, one per row
  MtxDbl initialIterates;

  /// best variables
  MtxDbl bestVars;

  /// best objective thus far
  double bestFunction;

};

} // end namespace nkm

#endif
