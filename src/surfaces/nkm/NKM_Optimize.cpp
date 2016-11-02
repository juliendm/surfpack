#include "NKM_Optimize.hpp"
#include "NKM_SurfPackModel.hpp"
#include <cfloat>
#include <cstdlib>

// define array limits hard-wired in DIRECT
// maxdim (same as maxor)
#define NCSU_DIRECT_MAXDIM 64
// maxfunc = 90000-20
#define NCSU_DIRECT_MAXFUNC 89980

#ifdef HAVE_CONFIG_H
// Tolerate F77_FUNC macro redefinition warnings in the autotools build
#define CONMIN_F77      F77_FUNC(conmin,CONMIN)
#define NCSU_DIRECT_F77 F77_FUNC_(ncsuopt_direct,NCSUOPT_DIRECT)

#else
// Use the CMake generated fortran name mangling macros (eliminate warnings)
#include "surf77_config.h"
#define CONMIN_F77      SURF77_GLOBAL(conmin,CONMIN)
#define NCSU_DIRECT_F77 SURF77_GLOBAL_(ncsuopt_direct,NCSUOPT_DIRECT)
#endif

extern "C" {

void CONMIN_F77(double* candidate, double* lowerb, double* upperb,
                double* constraint_values,
                double* scal, double* df, double* a, double* s, double* g1,
                double* g2, double* b, double* c,
                int* isc, int* ic, int* ms1,
                int& n1, int& n2, int& n3, int& n4, int& n5,
                double& delfun, double& dabfun, double& fdch, double& fdchm,
                double& ct, double& ctmin, double& ctl, double& ctlmin,
                double& alphax, double& abobj1, double& theta,
                double& obj,
                int& numdv, int& ncon, int& nside, int& iprint, int& nfdg,
                int& nscal, int& linobj, int& itmax, int& itrm, int& incdir,
                int& igoto, int& nac, int& info, int& infog, int& iter);

void NCSU_DIRECT_F77(int (*objfun)(int *n, double c[], double l[], double u[],
				   int point[], int *maxI, int *start,
				   int *maxfunc, double fvec[], int iidata[],
				   int *iisize, double ddata[], int *idsize,
				   char cdata[], int *icsize),
		     double* x, int& n, double& eps, int& maxf, int& maxT,
		     double& fmin, double* l, double* u, int& algmethod,
		     int& ierror, int& logfile, double& fglobal, double& fglper,
		     double& volper, double& sigmaper, int* idata, int& isize, 
		     double* ddata, int& dsize, char* cdata, int& csize,
		     int& quiet_flag);
}

namespace nkm {

OptimizationProblem* OptimizationProblem::optimizationProblemInstance(NULL);

// TODO: move to Teuchos, use putScalar (no need for bds check)

void OptimizationProblem::lower_bound(int i, double lb)
{ lowerBounds(i,0) = lb; }

void OptimizationProblem::upper_bound(int i, double ub)
{ upperBounds(i,0) = ub; }

void OptimizationProblem::initial_iterate(int i, double x0)
{ initialIterates(i,0) = x0; }

void OptimizationProblem::add_initial_iterates(MtxDbl& init_iterates_to_add)
{
  assert(init_iterates_to_add.getNRows()==numDesignVar);
  int numsofar=initialIterates.getNCols();
  int numtoadd=init_iterates_to_add.getNCols();
  //printf("numsofar=%d numtoadd=%d\n",numsofar,numtoadd);
  initialIterates.resize(numDesignVar,numsofar+numtoadd);
  MtxInt icols(numtoadd,1);
  for(int i=0; i<numtoadd; ++i)
    icols(i,0)=i+numsofar;
  initialIterates.putCols(init_iterates_to_add,icols);
}

// no treatment of constraints for now
void OptimizationProblem::multistart_conmin_optimize(int num_guesses)
{
  assert(num_guesses >= 1);

  MtxDbl guess(numDesignVar,1);
  double best_obj;
  bestFunction = DBL_MAX;

  double obj;
  MtxDbl con_out(10,1);

  // iterate over provided then possibly random guesses
  for (int iguess = 0; iguess < num_guesses; ++iguess) {
    //printf("iguess/num_guesses=%d/%d ",iguess,num_guesses);
    theModel.set_conmin_parameters(*this);

    retrieve_initial_iterate(iguess, guess);
    
    // TODO: put switch here for optimizer choice
    optimize_with_conmin(guess, best_obj);

    theModel.objectiveAndConstraints(obj,con_out,guess);
    //printf("} obj=%g\n",best_obj);


    // this will update guess to the best point found
    if(best_obj < bestFunction) {
      bestFunction = best_obj;
      bestVars = guess;
    }
  }

}


// single pass of CONMIN; if no initial iterate, use random guess
void OptimizationProblem::conmin_optimize()
{
  theModel.set_conmin_parameters(*this);
  // directly update bestVars/Functions in iteration
  retrieve_initial_iterate(0, bestVars);
  optimize_with_conmin(bestVars, bestFunction);
}


// DiRECT optimization within bounds, option of hidden condition
// number constraint
void OptimizationProblem::direct_optimize()
{
  theModel.set_direct_parameters(*this);
  // directly update bestVars/Functions in iteration
  optimize_with_direct(bestFunction);
}


//void OptimizationProblem::optimize_with_conmin(MtxDbl& guess)
void OptimizationProblem::optimize_with_conmin(MtxDbl& guess, 
					       double& final_val)
{
/*  The following was copied from the conmin user's manual found at 
    http://www.eng.buffalo.edu/Research/MODEL/mdo.test.orig/CONMIN/manual.html#Param_Main

SECTION IV
PARAMETERS DEFINED IN MAIN PROGRAM
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

NDV      <numDesignVar> Number of decision variables, X(I), contained in 
         vector X.

ITMAX    Default value = 10.  Maximum number of iterations in the
         minimization process.  If NFDG.EQ.0 each iteration requires one
         set of gradient computations (INFO = 3 or 4) and approximately
         three function evaluations (INFO = 1 or 2).  If NFDG.GT.0
         each iteration requires approximately NDV + 3 function
         evaluations (INFO = 1 or 2).

NCON     Number of constraint functions, G(J).  NCON may be zero.

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

NSCAL    Scaling control parameter.  The decision variables will be
         scaled linearly.

         NSCAL.LT.0:  Scale variables X(I) by dividing by SCAL(I), where
                      vector SCAL is defined by the user.

         NSCAL.EQ.0:  Do not scale the variables.

         NSCAL.GT.0:  Scale the variables every NSCAL iterations.
                      Variables are normalized so that scaled
                      X(I) = X(I)/ABS(X(I)).  When using this option, it
                      is desirable that NSCAL = ICNDIR if ICNDIR is input
                      as nonzero, and NSCAL = NDV + 1 in ICNDIR is input
                      as zero.

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

THETA    Default value = 1.0.  Not used if NCON = NSIDE = 0.  Mean value
         of the push-off factor in the method of feasible directions.
         A larger value of THETA is desirable if the constraints, G(J),
         are known to be highly nonlinear, and a smaller value may be
         used if all G(J) are known to be nearly linear.  The actual
         value of the push-off factor used in the program is a quadratic
         function of each G(J), varying from 0.0 for G(J) = CT to 4.0*THETA
         for G(J) = ABS(CT).  A value of THETA = 0.0 is used in the
         program for constraints which are identified by the user to be
         strictly linear.  THETA is called a "push-off" factor because
         it pushes the design away from the active constraints into the
         feasible region.  The default value is usually adequate.

PHI      Default value = 5.0.  Not used if NCON = NSIDE = 0.
         Participation coefficient, used if a design is infeasible
         (one or more G(J).GT.ABS(CT)).  PHI is a measure of how hard
         the design will be "pushed" towards the feasible region and
         is, in effect, a penalty parameter.  If in a given problem, a
         feasible solution cannot be obtained with the default value,
         PHI should be increased, and the problem run again.  If a
         feasible solution cannot be obtained with PHI = 100, it is
         probable that no feasible solution exists.  The default value
         is usually adequate.

NACMX1   Not used if NSIDE = NCON = 0.  1 plus user's best estimate of
         the maximum number of constraints (including side constraints,
         VLB(I) and VUB(I)) which will be active at any given time in
         the minimization process.  NACMX1 = number of rows in array A.
         If NAC + 1 ever exceeds this value, the minimization process will
         be terminated, an error message will be printed, and control
         will return to the main program.  NACMX1 will never exceed
         NDV + 1 if all constraints G(J) and bounds VLB(I) and VUB(I)
         are independent.  A reasonable value for NACMX1 (and the
         corresponding dimension of array A) is MIN(40, NDV + 1),
         where the minimum of 40 will only apply for large problems
         and is arbitrary, based on the observation that even for very
         large problems (over a hundred X(I) and several thousand G(J)),
         it is uncommon for many constraints to be active at any time
         in the minimization process (the optimum solution is seldom
         "fully constrained" for very large nonlinear problems).

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

LINOBJ   Not used if NCON = NSIDE = 0.  Linear objective function
         identifier.  If the objective, OBJ, is specifically known to
         be a strictly linear function of the decision variables, X(I),
         set LINOBJ = 1.  If OBJ is a general nonlinear function, set
         LINOBJ = 0.

ITRM     Default value = 3.  Number of consecutive iterations to indicate
         convergence by relative or absolute changes, DELFUN or DABFUN.

X(N1)    Vector of decision variables, X(I), I = 1, NDV.  The initial
         X-vector contains the user's best estimate of the set of optimum
         design variables.

VLB(N1)  Used only if NSIDE.NE.0.  VLB(I) is the lower allowable value
         (lower bound) of variable X(I).  If one or more variables, X(I),
         do not have lower bounds, the corresponding VLB(I) must be
          initialized to a very large negative number (say -1.0E+10).

VUB(N1)  Used only if NSIDE.NE.0.  VUB(I) is the maximum allowable value
         (upper bound) of X(I).  If one or more variables, X(I), do not
         have upper bounds, the corresponding VUB(I) must be initialized
         to a very large positive number (say 1.0E+10).

SCAL(N5) Not used if NSCAL = 0.  Vector of scaling parameters.  If
         NSCAL.GT.0 vector SCAL need not be initialized since SCAL will
         be defined in CONMIN and its associated routines.  If NSCAL.LT.0,
         vector SCAL is initialized in the main program, and the scaled
         variables X(I) = X(I)/SCAL(I).  Efficiency of the optimization
         process can sometimes be improved if the variables are either
         normalized or are scaled in such a way that the partial deri-
         vative of the objective function, OBJ, with respect to variable
         X(I) is of the same order of magnitude for all X(I).  SCAL(I)
         must be greater than zero because a negative value of SCAL(I)
         will result in a change of sign of X(I) and possibly yield
         erroneous optimization results.  The decision of if, and how, the
         variables should be scaled is highly problem dependent, and some
         experimentation is desirable for any given class of problems.

ISC(N8)  Not used if NCON = 0.  Linear constraint identification vector.
         If constraint G(J) is known to be a linear function of the
         decision variables, X(I), ISC(I) should be initialized to
         ISC(I) = 1.  If constraint G(J) is nonlinear ISC(I) is initialized
         to ISC(I) = 0.  Identification of linear constraints may improve
         efficiency of the optimization process and is therefore desirable,
         but is not essential.  If G(J) is not specifically known to be
         linear, set ISC(I) = 0.

SECTION V
PARAMETERS DEFINED IN EXTERNAL ROUTINE SUB1.

OBJ       Value of objective function for the current decision variables,
          X(I), I = 1, NDV contained in vector X.  Calculate OBJ if
          INFO = 1 or INFO = 2.

G(N2)     Not used if NCON = NSIDE = 0.  Vector containing all constraint
          functions, G(J), J = 1, NCON for current decision variables, X.
          Calculate G(J), J = 1, NCON if INFO = 2.

DF(N1)    Analytic gradient of the objective function for the current
          decision variables, X(I).  DF(I) contains the partial derivative
          of OBJ with respect to X(I).  Calculate DF(I), I = 1,
          NDV if INFO = 3 or INFO = 4 and if NFDG = 0 or NFDG = 2.

NAC       Number of active and violated constraints (G(J).GE.CT).
          Calculate NAC if INFO = 4 and NFDG = 0.

A(N4,N3)  Not used if NCON = NSIDE = 0.  Gradients of active or violated
          constraints, for current decision variables, X(I).
          A(J,I) contains the gradient of the Jth active or violated
          constraint, G(J), with respect to the Ith decision variable,
          X(I) for J = 1, NAC and I = 1, NDV.  Calculate if INFO = 4
          and NFDG = 0.  [SEE ADDENDUM]

IC(N4)    Identifies which constraints are active or violated.  IC(J)
          contains the number of the Jth active or violated constraint
          for J = 1, NAC.  For example, if G(10) is the first active
          or violated constraint (G(J).LT.CT, J = 1,9), set IC(1) = 10.
          Calculate if INFO = 4 and NFDG = 0.


CONMIN USER'S MANUAL ADDENDUM
by G. N. VANDERPLAATS
May, 1978
http://www.eng.buffalo.edu/Research/MODEL/mdo.test.orig/CONMIN/manual.html#Addendum

PROGRAM ORGANIZATION
The original version of CONMIN was written such that a user-supplied 
subroutine was called by CONMIN for function and gradient calculations. 
The current version is organized such that the function evaluation 
routine is contained in, or is called by, the main program. CONMIN 
executes according to the parameter IGOTO which must be initialized to 
zero. Figure A1 shows the required program organization. The purpose of 
this new logic is so that the program can be used in an overlay system 
or can be restarted in mid-execution.

ARRAY DIMENSIONS
The storage requirements for the arrays used in CONMIN have changed. The current storage requirements are:

DIMENSION X(N1), VLB(N1), VUB(N1), g(N2), SCAL(N1), df(N1), a(N1,N3), s(N1), g1(N2), g2(N2), b(N3,N3), c(N4), ISC(N2), ic(N3), ms1(N5)

where

    N1 = NDV + 2 
    N2 = NCON + 2*NDV 
    N3 = NACMX1 
    N4 = MAX (N3,NDV) 
    N5 = 2*N4
*/

  int N1 = numDesignVar + 2;
  int N2 = numConFunc + 2*numDesignVar;
  int N3 = numDesignVar + numConFunc + 1; //N3=NACMX1= 1 plus user's best estimate 
  //of the maximum number of constraints (including side constraints,
  //VLB(I) and VUB(I)) which will be active at any given time in
  //the minimization process. NACMX1 will never exceed
  //NDV + 1 if all constraints G(J) and bounds VLB(I) and VUB(I)
  //are independent.
  int N4 = N3;
  int N5 = 2*N4;
  
  int info = 0;                  ///CONMIN variable: status flag for optimization
  MtxDbl s(N1,1); s.zero();        ///Internal CONMIN array. Move direction in N-dimensional space.
  MtxDbl g1(N2,1); g1.zero();      ///Internal CONMIN array. Temporary storage of constraint values.  
  MtxDbl g2(N2,1); g2.zero();      ///Internal CONMIN array. Temporary storage of constraint values.
  MtxDbl B(N3,N3); B.zero();     ///Internal CONMIN array. Temporary storage of constraint values.
  MtxDbl c(N4,1);  c.zero();       ///Internal CONMIN array. Temporary storage for use with arrays B and S.
  MtxInt ms1(N5,1); ms1.zero();    ///Internal CONMIN array. Temporary storage for use with arrays B and S.
  MtxInt ic(N3,1); ic.zero();     ///Internal CONMIN array. Array of flags to identify active and violated constraints. I need to fill this in when I supply analytical gradients... see http://www.eng.buffalo.edu/Research/MODEL/mdo.test.orig/CONMIN/manual.html#List_4 for more details  
  double alphax = 0.1;           ///Internal CONMIN variable: 1-D search parameter.
  double abobj1 = 0.1;           ///Internal CONMIN variable: 1-D search parameter.
  double theta  = 1.0;           ///Internal CONMIN variable: mean value of push-off factor.
  //  double phi    = 5.0;           ///Internal CONMIN variable: "participation coefficient".
  int  nscal    = 0;             ///Internal CONMIN variable: scaling control parameter.
  MtxDbl scal(N1,1);  scal.zero(); ///Internal CONMIN array. Vector of scaling parameters for design parameter values.

  int  linobj   = 0;             ///Internal CONMIN variable: linear objective function identifier (unused).
  MtxInt isc(N2,1); isc.zero();    ///Internal CONMIN array. Array of flags to identify linear constraints. (not used in this implementation of CONMIN)
  
  int igoto = 0;                 ///Internal CONMIN variable: internal optimization termination flag. needs to be zero or CONMIN will carry over count of the number of evaluations of the objective function and its gradient from the construction of previous emulators, such as when the press metric of emulator quality is evaluated
  int nac   = 0;                 ///Internal CONMIN variable: number of active and violated constraints.  
  int infog = 0;                 ///Internal CONMIN variable: gradient information flag.
  int iter  = 0;                 ///Internal CONMIN variable: iteration count.

  ///conjugate direction restart parameter
  if(conminData.icndir==0) conminData.icndir=numDesignVar+1;

  MtxDbl query_pt(N1,1); //CONMIN CALLS THIS "X"
  MtxDbl lower_bounds(N1,1);
  MtxDbl upper_bounds(N1,1);
  theModel.makeGuessFeasible(guess,this); //need to find a better place to put this
  for(int ivar=0; ivar<numDesignVar; ivar++) {
    query_pt(ivar,0)=guess(ivar,0);
    lower_bounds(ivar,0)=lowerBounds(ivar,0);
    upper_bounds(ivar,0)=upperBounds(ivar,0);
  }

  double obj, dummy_obj;
  MtxDbl grad_obj(numDesignVar,1);
  MtxDbl con(numConFunc,1);
  MtxDbl grad_con(numConFunc,numDesignVar);

  MtxDbl df(N1,1); df.zero(); //the objective function gradient array that we need to pass into CONMIN (includes extra workspace), if a finite difference gradient is used this is a CONMIN internal array

  MtxDbl cv(N2,1); cv.zero(); //the constraint values array (called g(N2) in CONMIN's documentation) that we need to pass into CONMIN (includes extra workspace)

  MtxDbl A(N1,N3); A.zero(); //the gradients of constraints array that we need to pass into CONMIN (inludes extra workspace), if finite difference gradients are used this is a CONMIN internal array


  //assert((conminData.nfdg==0)||(conminData.nfdg==1)||(conminData.nfdg==2));
  int i, k;
  do {
    if(numConFunc>0) {
      //there are constraint FUNCTIONS
      //printf("  conmin iter=%d info=%d\n",iter,info);
      if(info>=2) {
	if(conminData.nfdg==1) {
	  //ConMin is requesting analyical GRADIENTS of the objective and constraint functions (but not the objective and constraint functions themselves)
	  theModel.objectiveAndConstraintsAndGradients(dummy_obj, con, 
						     grad_obj, grad_con, guess);
	  if(conminData.nfdg==1) {
	    nac=0;
	    for(k=0; k<numConFunc; k++) 
	      if(conminData.ct<=con(k,0)) {
		ic(nac,0)=k+1;
		for(i=0; i<numDesignVar; i++) 
		  A(i,nac) = grad_con(k,i);
		nac++;
	      }
	  }

	}
	else if(conminData.nfdg==2) //no analytical gradients for the constraints
	  theModel.objectiveAndGradient(dummy_obj, grad_obj, guess);
	
	for(i=0; i<numDesignVar; i++) 
	  df(i,0)=grad_obj(i,0);	
      }
      else{ //if(info==1) {
	//conmin is requesting the objective and constraint functions but NOT their gradients
	theModel.objectiveAndConstraints(obj, con, guess);
	for(k=0; k<numConFunc; k++) 
	  cv(k,0)=con(k,0);	  
      }
      //else
      //assert(info>0); //we shouldn't be able to get here
    }
    else{
      //there are NO constraint FUNCTIONS

      if((conminData.nfdg>0)&&(info>=2)) {
	//conmin is requesting the analytical GRADIENT of the objective function (but not the objective function itself)
	theModel.objectiveAndGradient(dummy_obj, grad_obj, guess);
	for(i=0; i<numDesignVar; i++) 
	  df(i,0)=grad_obj(i,0);
      }
      else{
	//conmin is requesting the objective function but NOT it's analytical gradient
	obj=theModel.objective(guess);
      }
    }

    CONMIN_F77(query_pt.ptr(0,0), lower_bounds.ptr(0,0), upper_bounds.ptr(0,0),
	       cv.ptr(0,0), scal.ptr(0,0), df.ptr(0,0), A.ptr(0,0), s.ptr(0,0),
	       g1.ptr(0,0), g2.ptr(0,0), B.ptr(0,0), c.ptr(0,0),
	       isc.ptr(0,0), ic.ptr(0,0), ms1.ptr(0,0), N1, N2, N3, N4, N5,
	       conminData.delfun, conminData.dabfun, 
	       conminData.fdch, conminData.fdchm,
	       conminData.ct, conminData.ctmin, conminData.ctl,
	       conminData.ctlmin, alphax, abobj1, theta, 
	       obj, numDesignVar, numConFunc, conminData.nside, 
	       conminData.iprint, conminData.nfdg, nscal, linobj, 
	       conminData.itmax, conminData.itrm, conminData.icndir, 
	       igoto, nac, info, infog, iter);

    for(i = 0; i<numDesignVar; i++) 
      guess(i,0) = query_pt(i,0);
    //printf("numDesignVar: %d query_pt.getNElems(): %d N1: %d\n",numDesignVar,query_pt.getNElems(),N1);
  } while (igoto != 0);

  final_val = obj;
}


// requires only the objective; TODO consider currVars as a guess
void OptimizationProblem::best_guess_optimize(int num_guesses)
{
  assert(num_guesses >= 1);

  MtxDbl guess(numDesignVar,1);
  bestFunction = DBL_MAX;

  // iterate over provided then possibly random guesses
  for (int iguess = 0; iguess < num_guesses; ++iguess) {
    retrieve_initial_iterate(iguess, guess);
    double obj = theModel.objective(guess);
    if(obj < bestFunction) {
      bestFunction = obj;
      bestVars = guess;
    }
  }

}


void OptimizationProblem::
optimize_with_direct(double& final_val)
{
  using std::cout;
  using std::cerr;

  if (directData.maxFunctionEvals > NCSU_DIRECT_MAXFUNC)
    std::cerr << "Error: Maximum function evaluations " 
	 << directData.maxFunctionEvals << "\nexceeds DiRECT algorithm limit " 
	 << NCSU_DIRECT_MAXFUNC << std::endl;
  if (numDesignVar > NCSU_DIRECT_MAXDIM)
    std::cerr << "Error: " << numDesignVar << " variables exceeds DiRECT algorithm "
	 << "limit of " << NCSU_DIRECT_MAXDIM << std::endl;
  if (directData.maxFunctionEvals > NCSU_DIRECT_MAXFUNC || 
      numDesignVar > NCSU_DIRECT_MAXDIM)
    std::exit(-1);

  // INITIALIZATION
  OptimizationProblem* prev_instance = optimizationProblemInstance;
  optimizationProblemInstance = this;

  int ierror, num_cv = numDesignVar, algmethod = 1, logfile = 13,
    quiet_flag  = directData.verboseOutput ? 0 : 1;
  double fmin = 0., eps = 1.e-4;

  // terminate when size of box  w/ f_min < sigmaper*size of orig box
  double sigmaper = 
    (directData.minBoxSize >= 0.) ? directData.minBoxSize : 1.e-4;
  // terminate when volume of box w/ f_min < volper*volume of orig box
  double volper = 
    (directData.volBoxSize >= 0.) ? directData.volBoxSize : 1.e-6;
  // convergence tolerance for target solution (DIRECT wants 0. when inactive)
  double fglper = 
    (directData.solutionTarget > -DBL_MAX) ? directData.convergenceTol : 0.;

  // for passing additional data to objective_eval()
  int isize = 0, dsize = 0, csize = 0;
  int*    idata = NULL;
  double* ddata = NULL;
  char*   cdata = NULL;

  NCSU_DIRECT_F77(direct_objective_eval, bestVars.ptr(0,0), num_cv, eps,
		  directData.maxFunctionEvals, directData.maxIterations, fmin, 
		  lowerBounds.ptr(0,0), upperBounds.ptr(0,0), algmethod, ierror, 
		  logfile, directData.solutionTarget, fglper, volper, sigmaper,
		  idata, isize, ddata, dsize, cdata, csize, quiet_flag);

  if (ierror < 0) {
    std::cerr << "NCSU DIRECT failed with fatal error code " << ierror << "\n";
    switch (ierror) {
    case -1:
      std::cerr << "(variable lower bounds must be strictly less than upper bounds)";
      break;
    case -2:
      std::cerr << "(maximum function evaluations is too large)";
      break;
    case -3:
      std::cerr << "(initialization in DIRpreprc failed)";
      break;
    case -4:
      std::cerr << "(error in creation of the sample points)";
      break;
    case -5:
      std::cerr << "(error occurred in sampling the function)";
      break;
    case -6:
      std::cerr << "(maximum iterations is too large)";
      break;
    default:
      std::cerr << "(unknown error code)";
    }
    std::cerr << "\nSee \"Calling DIRECT\" section in DIRECT Version 2.0 User Guide"
	 << ".\n" << std::endl;
    std::exit(-1);
  }
  else if (directData.verboseOutput) {
    std::cout << "NCSU DIRECT succeeded with code " << ierror << "\n";
    switch (ierror) {
    case 1:
      std::cout << "(maximum function evaluations exceeded)";
      break;;
    case 2:
      std::cout << "(maximum iterations reached)";
      break;;
    case 3:
      std::cout << "(prescribed global minimum reached within tolerance)";
      break;;
    case 4:
      std::cout << "(best rectangle reduced from original volume by prescribed "
	   << "fraction)";
      break;;
    case 5:
      std::cout << "(best rectangle measure is less than prescribed min box size)";
      break;;
    default:
      std::cout << "(unknown code)";
      break;;
    }
    std::cout << std::endl;
  }

  // FINALIZE
  optimizationProblemInstance = prev_instance;
  final_val = fmin;

  //std::cout << "fmin = " << fmin << "; theta = ";
  //for (int vi=0; vi<numDesignVar; ++vi)
  //std::cout << bestVars(vi) << " ";
  //std::cout << std::endl;

}


/// Modified batch evaluator that accepts multiple points and returns
/// corresponding vector of functions in fvec.  Must be used with modified
/// DIRECT src (DIRbatch.f).
int OptimizationProblem::
direct_objective_eval(int *n, double c[], double l[], double u[], int point[],
		      int *maxI, int *start, int *maxfunc, double fvec[],
		      int iidata[], int *iisize, double ddata[], int *idsize, 
		      char cdata[], int *icsize)
{
  int cnt = *start-1; // starting index into fvec
  int nx  = *n;       // dimension of design vector x.
  
  // number of trial points to evaluate
  // if initial point, we have a single point to evaluate
  int np = (*start == 1) ? 1 : *maxI*2;

  // loop over trial points, lift scaling, synchronously evaluate
  MtxDbl curr_vars(nx,1);
  int pos = *start-1; // only used for second eval and beyond
  for (int j=0; j<np; j++) {

    if (*start == 1)
      for (int i=0; i<nx; i++)
	curr_vars(1,0) = (c[i]+u[i])*l[i];
    else {
      for (int i=0; i<nx; i++) {
	// c[pos+i*maxfunc] = c(pos,i) in Fortran.
	double ci=c[pos+i*(*maxfunc)];
	curr_vars(i,0) = (ci + u[i])*l[i];
      }
      pos = point[pos]-1;
    }

    // choose between hidden constraint and unconstrained formula
    if (optimizationProblemInstance->directData.constraintsPresent) {
      
      double obj;
      MtxDbl con(optimizationProblemInstance->numConFunc,1);

      optimizationProblemInstance->
	theModel.objectiveAndConstraints(obj, con, curr_vars);

       // return function values
      fvec[cnt+j] = obj;

      // set flag to 1 if infeasible w.r.t. ANY constraint
      int infeasible = 0;
      //std::cout << "numConFunc=" << optimizationProblemInstance->numConFunc;
      for(int k=0; k<optimizationProblemInstance->numConFunc; k++) 
	if (!(con(k,0) < 0.0)) {
	  infeasible = 1;
	  //std::cout << "constraint violated" << std::endl;
	  break;
	}
      //if(infeasible==0)
      //std::cout << "a feasible solution exists" << std::endl;
      fvec[cnt+(*maxfunc)+j] = infeasible;

    }
    else {
      // return function values
      fvec[cnt+j] = optimizationProblemInstance->theModel.objective(curr_vars);
      // flag: successful eval
      fvec[cnt+(*maxfunc)+j] = 0; 
    }


  } // end evaluation loop over points

  return 0;
}


// return an initial iterate from the available list, otherwise return
// a random guess
void OptimizationProblem::retrieve_initial_iterate(int it_ind, MtxDbl& iterate)
{
  if (it_ind < initialIterates.getNCols()) {
    assert(initialIterates.getNRows() == numDesignVar);
    initialIterates.getCols(iterate, it_ind);
    //printf("opt:retrieve_initial_iterate it_ind=%d NCols=%d\n",it_ind,initialIterates.getNCols()); fflush(stdout);
  }
  else{
    //printf("opt:retrieve_initial_iterate getRandGuess\n"); fflush(stdout);
    getRandGuess(iterate);
  }
  //printf("iguess=%d EulAng=[ %g",iguess,guess(0,0));
  //for(int i=1; i<guess.getNElems(); i++)
  //printf(",%g ",guess(i,0));
  //printf("]\n");
}

void OptimizationProblem::getRandGuess(MtxDbl& guess) const
{
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  guess.newSize(numDesignVar,1);
  //printf("getRandGuess: lowerBounds.size=[%d %d] upperBounds=[%d %d]\n",
  //lowerBounds.getNRow(),lowerBounds.getNCols(),
  //upperBounds.getNRows(),upperBounds.getNCols()); fflush(stdout);
  for(int i=0;i<numDesignVar;i++)
    guess(i,0)=(std::rand() % mymod) *
      (upperBounds(i,0)-lowerBounds(i,0))/mymod + lowerBounds(i,0);
}

const MtxDbl& OptimizationProblem::best_point() const
{ return bestVars; }

}
