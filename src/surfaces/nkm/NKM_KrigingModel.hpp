/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_HPP__
#define __KRIGING_MODEL_HPP__

//#include "surfpack_system_headers.h"
#include "NKM_SurfPack.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_SurfPackModel.hpp"
#include "NKM_Optimize.hpp"
//#include "NKM_LinearRegressionModel.hpp"
#include <map>
#include <string>

namespace nkm {

typedef std::map< std::string, std::string> ParamMap;

// enumerated type stored in nkm::KrigingModel::corrFunc see below for more details
enum {DEFAULT_CORR_FUNC, GAUSSIAN_CORR_FUNC, EXP_CORR_FUNC, POW_EXP_CORR_FUNC, MATERN_CORR_FUNC};

// BMA TODO: Use more descriptive names for variables?

/** 
    KrigingModel: a class for creating and evaluating Gaussian process
    emulators with constant, linear, or quadratic trend function.
    Options for:

    * choice of optimizer
    * coordinate rotation
    * evaluation with gradients
    * nugget to control ill-conditioning.
    * optimal subset selection to control ill-conditioning
    * correlation Function (powered exponential or matern families)
*/
class KrigingModel: public SurfPackModel
{

public:

  // MtxDbl& makeGuessFeasible(MtxDbl& nat_log_corr_len, OptimizationProblem *opt);

  //returns a string indicating the correlation function
  std::string get_corr_func() const;

  std::string model_summary_string() const;

  // Return a string containing the model in "algebraic" format
  std::string asString() const;

  // BMA TODO: can we redesign so these need not be public?
  void set_conmin_parameters(OptimizationProblem& opt) const;

  void set_direct_parameters(OptimizationProblem& opt) const;

  // Creating KrigingModels

  /// Default constructor
  KrigingModel() : ifChooseNug(false), ifAssumeRcondZero(false), ifPrescribedNug(false), nug(0.0), XR(sdBuild.xr)
  { /* empty constructor */ };
  
  /// Standard KrigingModel constructor
  KrigingModel(const SurfData& sd, const ParamMap& params);

  /// After construction a Kriging model must be created with this
  /// function (TODO: add builtFlag for safety)
  void create();


  // Evaluating Kriging Models

  /// evaluate (y) the Kriging Model at a single point (xr is a Real row vector)
  double evaluate(const MtxDbl& xr);

  /// evaluate (y) the Kriging Model at a collection of points xr, one per row
  MtxDbl& evaluate(MtxDbl& y, const MtxDbl& xr);

  /// evaluate the KrigingModel's adjusted variance at a single point
  double eval_variance(const MtxDbl& xr);

  /// evaluate the KrigingModel's adjusted variance at a collection of points xr, one per row
  MtxDbl& eval_variance(MtxDbl& adj_var, const MtxDbl& xr);

  //double get_unadjusted_variance(){return (estVarianceMLE*scaler.unScaleFactorVarY());};
  
  /// evaluate the partial first derivatives with respect to xr of the models adjusted mean
  MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr);

  /// evaluate the partial second derivatives with respect to xr of the models adjusted mean... this gives you the lower triangular, including diagonal, part of the Hessian(s), with each evaluation point being a row in both xr (input) and d2y(output)
  MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr);

  // Helpers for solving correlation optimization problems

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "KrigingModel" good)  
  inline double objective(const MtxDbl& nat_log_corr_len) {
    MtxDbl corr_len(numTheta,1);
    for(int i=0; i<numTheta; ++i)
      corr_len(i,0)=std::exp(nat_log_corr_len(i,0));
    correlations.newSize(numTheta,1);
    get_theta_from_corr_len(correlations,corr_len);
    masterObjectiveAndConstraints(correlations, 1, 0);
    //printf("[objective]");
    return obj;
  };
    
  /// objective plus condition number constraints
  //void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
  inline void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
				      const MtxDbl& nat_log_corr_len) {
    //printf("entered objectiveAndConstraints\n");  fflush(stdout);
    MtxDbl corr_len(numTheta,1);
    for(int i=0; i<numTheta; ++i)
      corr_len(i,0)=std::exp(nat_log_corr_len(i,0));
    correlations.newSize(numTheta,1);
    get_theta_from_corr_len(correlations,corr_len);
    con_out.newSize(numConFunc,1);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i,0)=0.5*std::exp(-2.0*nat_log_corr_len(i,0));
    //printf("about to enter masterObjectiveAndConstraints\n"); fflush(stdout);
    masterObjectiveAndConstraints(correlations, 1, 1);
    //printf("left masterObjectiveAndConstraints\n"); fflush(stdout);
    obj_out=obj;
    for(int i=0; i<numConFunc; i++){
      //printf("i=%d ",i); fflush(stdout);
      con_out(i,0)=con(i,0);
    }
    //con_out.copy(con);
    //printf("[objectiveAndConstraints]");
    //printf("leaving objectiveAndConstraints\n");  fflush(stdout);
    return;
  };

  /// return the Number of Trend functions, the trend is represented by an
  /// arbitrary order multidimensional polynomial, individual trend functions
  /// are the separate additive terms in that multidimensional polynomial
  inline int getNTrend() const
  { return (Poly.getNCols());   } 

  // return the likelihood of this model
  inline double getLikelihood()
  { return likelihood; }

  static int min_coefficients(int nvars, int poly_order) 
  {
    return num_multi_dim_poly_coef(nvars,poly_order)+nvars;
  };

  void getRandGuess(MtxDbl& guess) const;

private:
  
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class SurfPoint data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

  // helper functions
  void preAllocateMaxMemory();
  void reorderCopyRtoRChol();
  void nuggetSelectingCholR();
  void equationSelectingCholR();
  void trendSelectingPivotedCholesky();

  /// this function calculates the objective function (negative log
  /// likelihood) and/or the constraint functions and/or their analytical
  /// gradients and/or the hessian of the objective function using a 
  /// precompute and store (store across sequential calls to this function) 
  /// strategy to reduce the computational cost.  To ensure that precomputed
  /// and stored values are not changed externally this function return no
  /// output, instead member variables must be copied out by wrapper functions.
  /// The objective and contraint derivative modes are bit flags, i.e. 
  /// each is the sum of 2^(all orders of derivative you want). KRD 2010.05.13
  void masterObjectiveAndConstraints(const MtxDbl& theta, int obj_der_mode, 
				     int con_der_mode);

  //void set_conmin_parameters(OptimizationProblem& opt) const;

  /// evaluate the trend function g(xr), using class member Poly
  inline MtxDbl& eval_trend_fn(MtxDbl& g, const MtxDbl& xr) {
    return (evaluate_poly_basis(g, flyPoly, Poly, xr));
  }

  inline MtxDbl& eval_der_trend_fn(MtxDbl& dg, const MtxInt& der, 
				   const MtxDbl& xr) {
    return (evaluate_poly_der_basis(dg, flyPoly, derivBetaHat, Poly, der, xr));
    return dg;
  }


  //the following matern_1pt5_... and matern_2pt5_... functions don't really need to be member functions as they don't access any data members

  /// multiply exponential corr func by this to get matern 1.5 corr func
  inline double matern_1pt5_coef(double theta_abs_dx) const {
    return 1.0+theta_abs_dx;
  };
  /// multiply matern 1.5 corr func by this to get d1 of matern 1.5 corr func
  inline double matern_1pt5_d1_mult_r(double theta, double dx) const {
    return -theta*theta*dx/matern_1pt5_coef(theta*std::fabs(dx));
  };
  /** multiply matern 1.5 corr func by this to get d2 of matern 1.5 corr func
      1D MATERN_CORR_FUNC 1.5 r(x1,x2) is twice+ differential except 
      where x1==x2 this is correct for x1!=x2 */
  inline double matern_1pt5_d2_mult_r(double theta, double dx) const{
    return theta*theta*(1.0-2.0/matern_1pt5_coef(theta*std::fabs(dx)));
  };


  /// multiply exponential corr func by this to get matern 2.5 corr func
  inline double matern_2pt5_coef(double theta_abs_dx) const{
    return 1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx/3.0;
  };
  /// multiply matern 2.5 corr func by this to get d1 of matern 2.5 corr func
  inline double matern_2pt5_d1_mult_r(double theta, double dx) const{
    double theta_abs_dx=theta*std::fabs(dx);
    return (-theta*theta*dx*(1.0+theta_abs_dx)/
	    (3.0*matern_2pt5_coef(theta_abs_dx)));
  };
  /// multiply matern 2.5 corr func by this to get d2 of matern 2.5 corr func
  inline double matern_2pt5_d2_mult_r(double theta, double dx) const{
    double theta_abs_dx=theta*std::fabs(dx);
    return -theta*theta*(1.0+theta_abs_dx-theta_abs_dx*theta_abs_dx)/
      (3.0*matern_2pt5_coef(theta_abs_dx));
  };


  // BMA TODO: these docs need updating

  /** converts from correlation lengths to theta
      for powered exponential (including exponential and Gaussian)
          theta=1/(powExpCorrLenPow*corr_len^powExpCorrLenPow)
      for matern (excluding Gaussian)
          theta= sqrt(2*maternCorrFuncNu)/corr_len  */
  MtxDbl& get_theta_from_corr_len(MtxDbl& theta, const MtxDbl& corr_len) const;
  /** converts from theta to correlation lengths 
      for powered exponential (including exponential and Gaussian)
          theta=1/(powExpCorrLenPow*corr_len^powExpCorrLenPow)
      for matern (excluding Gaussian)
          theta= sqrt(2*maternCorrFuncNu)/corr_len  */
  MtxDbl& get_corr_len_from_theta(MtxDbl& corr_len, const MtxDbl& theta) const;

  MtxDbl& eval_kriging_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const;
  MtxDbl& eval_gek_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const;
  /** r(i,j)=corr_func(xr(i,:),XR(j,:);theta(:)) choices for correlation 
      function are gaussian, exponential, powered exponential with 1<power<2, 
      and matern with nu=1.5 or 2.5 (gaussian and exponential are pulled out
      for efficient implementation and because they belong to both families).
      (note that only matern 1.5, matern 2.5, and gaussian are avaliable for 
      Gradient Enhanced Kriging) The convention is that capital matrices are 
      for the data the model is built from, lower case matrices are for 
      arbitrary points to evaluate the model at. This function calls either
      the eval_kriging_* or eval_gek_* version of the same function depending 
      on wheter Kriging or Gradient Enhanced Kriging (GEK) is being used) */
  inline MtxDbl& correlation_matrix(MtxDbl& r, const MtxDbl& xr) const {
    if(buildDerOrder==0)
      return eval_kriging_correlation_matrix(r,xr);
    else if(buildDerOrder==1)
      return eval_gek_correlation_matrix(r,xr);
    else{
      std::cerr << "unsupported derivative order in\n  inline MtxDbl& correlation_matrix(MtxDbl& r, const MtxDble& xr) const\n";
      assert(false);
    }
  };

  MtxDbl& eval_kriging_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const;
  MtxDbl& eval_gek_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const;
  /** if r(i,j) is the evaluation of the corr_func(XR(:,i),xr(:,j)) then this
      function returns the matrix dr which is the derivative of matrix r with
      respect to dimension Ider of xr i.e.
      dr(i,j)=d(r(i,j))/d(xr(Ider,j)) 
      combining repeated calls can be used to get arbitrary (mixed) higher 
      order derivatives BUT this doesn't work when the two or more 
      derivatives are with respect to the same input variable; any unmixed 
      higher order derivatives (including when it is a component of an even 
      higher order mixed order derivative) requires special treatment. This 
      function calls eitherthe eval_kriging_* or eval_gek_* version of the 
      same function depending on wheter Kriging or Gradient Enhanced Kriging 
      (GEK) is being used) */
  inline MtxDbl& dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
					 const MtxDbl& xr, int Ider) const
  {
    if(buildDerOrder==0)
      return eval_kriging_dcorrelation_matrix_dxI(dr, r, xr, Ider);
    else if(buildDerOrder==1)
      return eval_gek_dcorrelation_matrix_dxI(dr, r, xr, Ider);
    else{
      std::cerr << "unsupported derivative order in\n inline MtxDbl& dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const\n";
      assert(false);
    }
  };
  
  MtxDbl& eval_kriging_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const;
  MtxDbl& eval_gek_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const;
  /** d2r(i,j)= d^2r(i,j)/dxr(Ider,j)dxr(Jder,j) where j is the point 
      index for xr, i is the point index for XR, and Ider and Jder are 
      the dimensions with respect to which the derivatives are being 
      taken, drI is the first derivative with respect to Ider (which 
      must be precomputed and passed  in), this function calls either 
      the eval_kriging_* or eval_gek_* version of the same function 
      depending on whether Kriging OR Gradient Enhanced Kriging (GEK) 
      is being used. */
  inline MtxDbl& d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const
  {
    if(buildDerOrder==0)
      return eval_kriging_d2correlation_matrix_dxIdxJ(d2r,drI,r,xr,Ider,Jder);
    else if(buildDerOrder==1)
      return eval_gek_d2correlation_matrix_dxIdxJ(d2r,drI,r,xr,Ider,Jder);
    else{
      std::cerr << "unsupported derivative order in\ninline MtxDbl& d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const\n";
      assert(false);
    }
  };

  /** R(i,j)=corr_func(XR(i,:),XR(j,:);theta(:)) where choices for for 
      correlation function are gaussian, exponential, powered exponential 
      with 1<power<2, and matern with nu=1.5 or 2.5 (gaussian and exponential 
      are pulled out for efficient implementation and because they belong 
      to both families).  All of the preceeding correlation functions are 
      available for Kriging, only Gaussian, Matern 1.5 and Matern 2.5 are
      available for GEK.  All correlation functions are implemented as
      Kriging R=something.*exp(Z*theta) (with reshapes) the something 
      depends on the correlation function (for gaussian, exponential, and 
      powered exponential that something is 1, and the multiplication by 
      1 is not actually done for the sake of efficiency) for the matern 
      function that something is matern_1pt5_coef or matern_2pt5_coef) of 
      course the definition of Z and theta differs for different correlation 
      functions.  Note that Z only stores what is needed to compute the 
      strictly lower (BELOW the diagonal) part of the Kriging R to save 
      memory and computation.  The Kriging R is symmetric and has ones on 
      the diagonal. The GEK R matrix is blocked into (1+numVarsr) by
      (1+numVarsr) submatrices.  Each submatric has numPoints by numPoints
      elements.  The upper-left-most submatrix is the Kriging R matrix
      the others submatrices are derivatives of the Kriging R with respect
      to the various dimensions of the first (rows) and second (columns) 
      inputs of the correlation function */
  void correlation_matrix(const MtxDbl& corr_vec);

  /** this function applies the nugget to the R matrix (a member variable)
      and stores the result in R (another member variable), i.e. it adds 
      nug to the diagonal of R. The convention is that capital matrices 
      are for the data the model is built from, lower case matrices are for 
      arbitrary points to evaluate the model at.  Once the nugget is added
      R is not strictly a correlation matrix */
  void apply_nugget_build();

  /** the Z matrix, Z=Z(XR), its definitition depends on the correlation 
      function
          for the gaussian correlation function 
              Z(ij,k)=-(XR(i,k)-XR(j,k))^2, 
          for the exponetial and matern correlation functions
              Z(ij,k)=-|XR(i,k)-XR(j,k)|     
	  for the powered exponential correlation function
	      Z(ii,k)=-|XR(i,k)-XR(j,k)|^powExpCorrFuncPow
	      where 1<powExpCorrFuncPow<2
      the Z matrix facilitates the efficient evaluation of the correlation 
      matrix R... R=something.*exp(Z*theta) 
      note that Z only holds what is needed to compute the strictly lower
      (below the diagonal) portion of R to save memory and computation.
      R is symmetric and has ones on the diagonal.  The convention is that 
      capital matrices are for the data the model is built from, lower 
      case matrices are for arbitrary points to evaluate the model at, 
      the Z and XR matrices are member variables so they don't need to be 
      passed in */
  MtxDbl& gen_Z_matrix();
  
  /** the order of the derivatives this Kriging Model was built for
      buildDerOrder=0  means function values only (regular Kriging)
      buildDerOrder=1  means function values plus gradients (Gradient Enhanced 
                       Kriging)
      buildDerOrder>=2 is currently not allowed, a later developer could 
                       implement Hessian Enhanced Kriging but KRD did not 
		       do this  */
  short buildDerOrder; 
  
  /** number of derivatives used to construct the Kriging/GEK model
      the zeroth-order derivative, i.e. function value itself counts
      as a derivative, so for Kriging nDer=1, for GEK nder=1+numVarsr */
  int nDer;

  /** a "Poly" style matrix (see "Poly" below) that stores derivatives 
      orders used to construct the Kriging/GEK model.  Der is a 
      numVarsr by nDer matrix. For Kriging Der is a numVarsr by 1 matrix 
      of zeros.  For GEK Der is a numVarsr by 1+numVarsr matrix whose 
      first (index zero) column is all zeros and columns with index 1 
      through numVarsr hold the identiy matrix which is the mixed partial
      derivative order representation of the gradient. */
  MtxInt Der;

  /** stores which correlation function we using, major choices are 
          powered exponential with 1<=power<=2 and 
          matern with nu=0.5,1.5,2.5 or "infinity" 
      There are 2 special cases that belong to both families and are 
      pulled out for efficient implementation these are
          gaussian correlation function
              equals powered exponential with power=2 
              equals matern with nu="infinity"
          exponential correlation function 
              equals powered exponential with power=1
	      equals matern with nu = 0.5
      after these are pulled out we have 
      powered exponential with 1<power<2 and
      matern with nu = 1.5 or 2.5
      corrFunc stores an enumerated type
  */
  short corrFunc;

  ///the power of the powered exponential family of correlation functions 
  double powExpCorrFuncPow;
  
  /// the "nu" parameter of the Matern family of correlation functions
  double maternCorrFuncNu;

  /** used to determine the "small feasible region" that we need to search to
      find good correlation lengths for the chosen correlation function */
  double aveDistBetweenPts;

  /// the upper bound of the small feasible region of correlation lengths
  double maxNatLogCorrLen;

  /// the lower bound of the small feasible region of correlation lengths
  double minNatLogCorrLen;

  /** the chosen natural log of the correlation LENGTHS (NOT correlation 
      PARAMETERS) */
  MtxDbl natLogCorrLen;

  /** the vector of correlation parameters (NOT correlation LENGTHS), these
      are the values determinened by the maximum likelihood optimization, the
      temporary in process version is called "theta" */
  MtxDbl correlations;

  /// NUMber of Real input VARiableS => NUMRVARS => NUMVARSR => numVarsr
  int numVarsr;

  /// NUMber of THETA for Real input variables... this is the number of correlation parameters (for real input variables), for Kriging numTheta=numVarsr, for radial basis functions numTheta=1 (radial basis functions, intended to be a derived class, are independent of direction)
  int numTheta;

  /** what optimization method did the user want us to use to determine the 
      correlation lengths of the Gaussian Process error model */
  std::string optimizationMethod;

  /** did the user specify correlation Lengths (either to use directly or 
      to use as the starting location of local optimization) */
  bool ifUserSpecifiedCorrLengths;

  /// number of starting locations for (possibly multistart) local optimization to try
  int numStarts;

  /// maximum number of sets of roughness parameters to try
  int maxTrials;

  ///used if optimization_method = global_local
  int maxTrialsGlobal; 

  ///used if optimization_method = global_local
  int maxTrialsLocal; 

  //"rcond" is now the only allowed constraint type, the eig option was removed
  //std::string constraintType;

  /** the number of constraint FUNCTIONS (typically these are nonlinear), 
      this number does NOT include box edge constraints for the inputs,
      those are handled separately, the method of computing analytical 
      derivatives of eigenvalues was questionable also the eigenvalue 
      approach was for the 2 norm condition number of the R matrix when 
      what we care about is actually the 1 norm condition number of the 
      R matrix, so the option for analytical constraints were removed which 
      means numConFunc should be exactly 1 */
  int numConFunc; 

  /** the correlation matrix R is considered to be "ill-conditioned" if
      if it's condition number exceeds this value. The constraint is that
      R not be ill conditioned */
  double maxCondNum;

  /** ifChooseNug==true tells KrigingModel to choose the smallest nugget it 
      needs to fix ill conditioning, ifChooseNug=0 tells KrigingModel not 
      to choose one (the default value for the nugget is zero) but the 
      user still has the option to prescribe a nugget the Nugget must be 
      a positive number */
  bool ifChooseNug; 
  
  /** iff ifChooseNug==true then ifAssumeRcondZero matters
      * for normal operations (ifAssumeRcondZero==false) a LAPACK Cholesky 
        factorization is performed, from that rcondR is calculated
      * if ifAssumeRcondZero==true then we skip this first LAPACK Cholesky and 
        assume rcondR=0.0 (to speed things up)
      If rcondR says R is ill conditioned, then the minimum size nugget that 
      is guaranteed to fix worst case ill-conditioning for the rcondR is 
      chosen, that nugget is applied to the diagonal of R, and a "second" 
      LAPACK Cholesky (this time of R with the nugget) is performed.
      ifAssumeRcondZero==true is a way to speed things up by always adding 
      a still very small nugget, it can be particulary useful for Gradient
      Enhanced Kriging if you would like to add a nugget.*/
  bool ifAssumeRcondZero;

  /** if ifPrescribedNug==true then the user has prescribed a nugget, 
      (think of the nugget a measurement noise term, it should be 
      roughly the variance of the measurement noise divided by the 
      variance of the output at the build data points) this will 
      cause the GP to smooth the data.  Typically this will be more 
      than sufficient to handle ill conditioning, but if not the bound 
      on rcond will be used to restrict the allowable range of correlation
      parameters so that all data points will be used */
  bool ifPrescribedNug;

  /** the nugget value sets the ammount of smoothing (approximation instead
      of interpolation) that the KrigingModel will use, it can also be used
      to fix ill conditioning, setting ifChooseNug=1 tells KrigingModel to 
      choose the smallest nugget needed to fix ill conditioning */
  double nug;

  /// the number of build points available
  int numPoints;

  /** which points are we keeping, for when using Pivoted Cholesky to select 
      an optimal subset of points to retain, if we're not selecting a subset 
      of points then this is ALL points.  When GEK is used all but the last
      point is required to be a whole point */
  MtxInt iPtsKeep; 

  /** the number of points retained after using Pivoted Cholesky to select
      an optimal subset of points. For GEK partial points are included in 
      this number */
  int numPointsKeep;
 
  /** only meaningful is GEK is used, 
      for Kriging numWholePointsKeep=numPointsKeep
      for GEK numWholePointsKeep is either numPointKeep or numPointsKeep-1 */
  int numWholePointsKeep; 

  /** only meaningful if GEK is used, the number of derivatives the last 
      retained point */
  int numExtraDerKeep; 

  /** the number of equations available to build the Kriging Model
      for regular Kriging numEqnAvail=numPoints;
      for Gradient Enhanced Kriging (GEK) numEqnAvail=(1+numVarsr)*numPoints */
  int numEqnAvail;

  /** the number of rows (and columns) in the R matrix.  For Kriging this
      is identical to numPointsKeep.  For GEK this is 
      (1+numVarsr)*numWholePointsKeep +
      (numPointsKeep-numWholePointsKeep)*(1+numExtraDerKeep) */
  int numRowsR;

  /** do we have an anchor point, i.e. one point that the user has required us
      to retain when we select a subset of points to build our Kriging model 
      from */
  bool ifHaveAnchorPoint;

  /// if we have an Anchor point what is its index?
  int  iAnchorPoint;

  /** the input the model was constructed from; convention is capital
      matrices are data model is built from, lower case matrices are
      arbitrary points to evaluate model at, using XR instead of X in
      anticipation of mixed real and integer input, use a reference XR 
      as shortcut/shorthand-notation to xr in the surfdata
  */
  MtxDbl& XR; 
  MtxDbl XRreorder;  //a reordered (by pivoted cholesky subset selection)
  //version of XR to make emulator EVALUATION fast

  /** the output at ALL available build data points, reshaped to a vector. 
      If GEK is used the function value and derivatives at a point are 
      sequential (i.e. a "whole point" at a time) */
  MtxDbl Yall;

  /** the likely reorderd subset of build point output data that the model 
      was constructed from. If GEK is used, all but the last point is 
      guaranteed to be a whole point (a function value immediately followed
      by the entire gradient in standard gradient order); the last retained
      point can be partial (i.e. can be missing some or all derivatives, but is 
      guaranteed to have the function value, the derivatives in the final 
      gradient were not reorderd before trailing entries were dropped). The
      convention is capital matrices are data model is built from, lower case
      matrices are arbitrary points to evaluate model at */
  MtxDbl Y;   

  /** the trend basis functions evaluated at all available build data points
      in their original order.  It has npoly rows. For Kriging it has numPoints
      columns.  For GEK it has (1+numVarsr)*numPoints columns with each "whole
      point" appearing as as 1+numVarsr sequential columns */
  MtxDbl Gall;

  /** the transpose of the matrix of trend function evaluations at the 
      (likely reorderd) subset of points used to build the Kriging Model.  
      If GEK is used, it will also contain derivatives of the original 
      polynomial basis functions. The convention is that capital matrices 
      are for the data the model is built from, lower case matrices are 
      for arbitrary points to evaluate the model at */
  MtxDbl Gtran;

  /** true if the user said to use only the main effects (no interaction/mixed 
      terms) in the polynomial basis */
  bool ifReducedPoly; 

  /** this is what the user asked for, highest total order of any term in the 
      polynomial basis */
  int polyOrderRequested; 
  
  /** this is what was actually used, can be less than what the user asked for 
      if the correlation matrix was ill-conditioned and we had to drop points 
      to fix the ill-conditioning and had to drop trend order to make it less 
      than the remaining points. */
  int polyOrder; 

  /** the number of equations needed for trend functions of order 0 through 
      polyOrderRequested */
  MtxInt numTrend; 

  /// the number of terms in the trend function that was actually used 
  int nTrend; 

  /** the indexes of the subset trend basis functions that a pivoted 
      Cholesky factorization of G*R^-1*G^T selected for retention, these
      indexes must be in "logical order" (monotonically increasing)*/
  MtxInt iTrendKeep;

  /** the polynomial powers for individual dimensions in each (trend) "basis 
      function" (for now the only choice of polynomial basis functions are
      multidimensional monomials) in a multidimensional polynomial of 
      arbitrary order
      For example if a multidimensional monomial is 
          coef * x0^0 * x1^2 * x2^0 * x3^1 *x4^0
      It's "Poly" representation would be 
	  [0 2 0 1 0]^T (matlab notation) with an associated coefficient 
	  stored in betaHat
      Each monomial in a polynomial is a column of the "Poly" matrix
      This format has many features that make it good for a "permanent"
      record of polynomials, for example it is convenient for taking 
      analytical derivatives.  But it is not particularly efficient to 
      evaluate when a set of numPoints points is stored as a matrix with 
      numVarsr rows and numPoints columns (where each column is a point).
  */
  MtxInt Poly;  

  /** flyPoly is "work space" for an on the fly generated 
      represenation of polynomials that is concise/fast to evaluate,
      particularly when a set of numPoints points is stored as a
      matrix with numVarsr rows and numPoints columns (where each column
      is a point).  But the "flypoly" style representation is not a
      convenient format for performing analytical derivatives of polynomials.
      However, converting from a "poly" representation to a "flypoly" 
      representation is trivially easy. See comments in NKM_SurfPack.hpp 
      for a description of the "flypoly" format. */
  MtxInt flyPoly; 

  /** the vector of coefficients of the trend functions (unadjusted mean)
      betaHat=(G*R^-1*G^T)^-1*(G*R^-1*Y) i.e. the generalized by R^-1 
      least squares fit, the generalization makes it unbiased */
  MtxDbl betaHat;

  /// modified coefficients for on the fly derivatives of the trend functions
  MtxDbl derivBetaHat;

  /** the Z matrix, Z=Z(XR), 
      * for the Gaussian Correlation function
        Z(ij,k)=-(XR(i,k)-XR(j,k))^2
      * for the Exponential, Matern 1.5 and Matern 2.5 Correlation functions
        Z(ij,k)=-|XR(i,k)-XR(j,k)|
      * for the powered Exponential Correlation function (other than the 
        Exponential and Gaussian correlation functions)
        Z(ij,k)=-(|XR(i,k)-XR(j,k)|^powExpCorrFuncPow)
      The Z matrix facilitates the efficient evaluation of the correlation 
      matrix R and its derivatives with respect to theta (the vector of 
      correlation parameters);
      i.e. R=coefficient*exp(Z*theta) the coefficient is 1.0 for the powered
      Exponential family of correlation functions, it is not 1.0 for the 
      Matern 1.5 and Matern 2.5 correlation functions.  The convention is 
      that capital matrices are for the data the model is built from, 
      lower case matrices are for arbitrary points to evaluate model at 
      size(Z)=[numVarsr nchoosek(numPoints,2)] */
  MtxDbl Z; 

  /** working memory (so we don't constantly need to allocate and deallocate it)
      used during the calculation of R, equals Z^T*theta (matrix multiplication 
      is used), 
      size(Ztran_theta)=[nchoosek(numPoints,2) 1] */
  MtxDbl Ztran_theta;

  /** deltaXR(ij,k)=XR(i,k)-XR(j,k) for the strictly lower triangular part of R
      this is useful for efficient computation of derivative enhanced R matrix 
      for Gradient Enhanced Kriging, if Kriging is used this will be an
      empty matrix (i.e. space will not be allocated for it) 
      size(deltaXR)=[nchoosek(numPoints,2) numVarsr] 
      having a transposed order of Z should make it faster to compute the GEK R
      matrix sine it is "blocked" into (1+numVarsr) by (1+numVarsr) submatrices.
      The size of each submatrix is numPoints by numPoints (i.e. the size of 
      the Kriging R matrix) */
  MtxDbl deltaXR;

  /** the "correlation matrix," for either regular Kriging or Gradient Enhanced
      Kriging, after possible inclusion of a nugget, use of a nugget causes 
      the KrigingModel to smooth i.e. approximate (which is useful if you would
      like to account for known measurement noise) rather than interpolate and 
      can be used to fix ill-conditioning. Technically R is only a Correlation
      Matrix if all of it's diagonal elements = 1.0, this is not the case if
      you've added a nugget or if you're using Gradient Enhanced Kriging */
  MtxDbl R;

  //we will need Rinv if we want to evaluate the INTEGRAL of the adjusted variance; if that gets implemented we should compute Rinv as the "last step" of the construction of the Kriging/GEK model
  //MtxDbl Rinv; 
  //Rinv_Gtran*inv(G_Rinv_Gtran)*Rinv_Gtran^T would also be needed for analtyical integration of the adjusted variance, we should calculate this ONCE after the optimization in complete (when we clear stuff, or maybe we should just go ahead and calculate the integral (mean) of adjusted variance and store the answer in case anyone ever asks for it, then we wouldn't have to store Rinv or the longer matrix forever)

  /** The lower triangular part of the Cholesky decomposition of R (the 
      correlation matrix after possible modification by the inclusion of 
      a nugget).  Keep this around to evaluate the adjusted variance. 
      The convention is that capital matrices are for the data the model 
      is built from, lower case matrices are for arbitrary points to 
      evaluate the model at */
  MtxDbl RChol;

  /** working memory for the factorization of R so that the equilibrated
      Cholesky Lapack wrapper won't have to allocate memory each time it is
      called, done for computational efficiency */
  MtxDbl scaleRChol;

  /** working memory for efficient computation of rcond when using pivoted 
      Choleksy to select an optimal subset of points to retain */
  MtxDbl sumAbsColR;

  /** working memory for efficient computation of rcond when using pivoted 
      Choleksy to select an optimal subset of points to retain */
  MtxDbl oneNormR;

  /** working memory for the bisection search for the maximum number of 
      points that can be retained after using Pivoted Cholesky to peform
      an optimal (in terms of unique information content) reordering 
      of points/equations, used when using Pivoted Cholesky to select an 
      optimal subset of points to retain */
  MtxDbl lapackRcondR;

  /** working memory for efficient computation of rcond when using pivoted 
      Choleksy to select an optimal subset of points to retain */
  MtxDbl rcondDblWork;

  /** working memory for efficient computation of rcond when using pivoted 
      Choleksy to select an optimal subset of points to retain */
  MtxInt rcondIntWork;

  /** the rcond (estimated reciprocal of the condition number) of the 
      modified correlation matrix, R */
  double rcondR; 

  /// rcond of (G^T*R^-1*G)
  double rcond_G_Rinv_Gtran;

  /// keep around to evaluate adjusted variance
  MtxDbl Rinv_Gtran;

  /// don't keep arround
  MtxDbl G_Rinv_Gtran;

  /// keep around to evaluate adjusted variance
  MtxDbl G_Rinv_Gtran_Chol;

  /** working memory used to calculate G_Rinv_Gtran_Chol efficiently (so we 
      don't have to constantly allocate/deallocate memory) */
  MtxDbl G_Rinv_Gtran_Chol_Scale;

  /** working memory used to calculate G_Rinv_Gtran_Chol efficiently (so we 
      don't have to constantly allocate/deallocate memory) */
  MtxDbl G_Rinv_Gtran_Chol_DblWork;

  /** working memory used to calculate G_Rinv_Gtran_Chol efficiently (so we 
      don't have to constantly allocate/deallocate memory) */
  MtxInt G_Rinv_Gtran_Chol_IntWork;

  /// working memory G*R^-1*Y
  MtxDbl G_Rinv_Y;

  /// working memory eps=epsilon=Y-G^T*betaHat
  MtxDbl eps;

  /** rhs = right hand side, rhs=Rinv*(Y-G^T*betaHat); The convention is 
      that capital matrices are for the data the model is built from, lower
      case matrices are for arbitrary points to evaluate the model at */
  MtxDbl rhs; 

  /// need keep around to evaluate adjusted variance
  double estVarianceMLE;

  /// the per equation log(likelihood) of the Kriging Model
  double likelihood; 

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int prevObjDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int prevConDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  MtxDbl prevTheta; //(numTheta,1)

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int maxObjDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int maxConDerMode;

  /// the objective function for the optimization of correlation lengths, it's the negative "per equation" log likelihood function
  double obj;

  /// the vector (a numConFunc by 1 matrix) of constraint functions for the optimization of the correlation lengths, it only needs to be a vector for compatibility with the nkm::OptimizationProblem class, otherwise it could be a double
  MtxDbl con; 
};

} // end namespace nkm

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template< class Archive >
void nkm::KrigingModel::serialize(Archive & archive, 
			  const unsigned int version)
{  

  archive & boost::serialization::base_object<nkm::SurfPackModel>(*this);
  archive & buildDerOrder;
  archive & nDer;
  archive & Der;
  archive & corrFunc;
  archive & powExpCorrFuncPow;
  archive & maternCorrFuncNu;
  archive & aveDistBetweenPts;
  archive & maxNatLogCorrLen;
  archive & minNatLogCorrLen;
  archive & natLogCorrLen;
  archive & correlations;
  archive & numVarsr;
  archive & numTheta;
  archive & optimizationMethod;
  archive & ifUserSpecifiedCorrLengths;
  archive & numStarts;      
  archive & maxTrials;
  archive & maxTrialsGlobal;
  archive & maxTrialsLocal;
  archive & numConFunc;
  archive & maxCondNum;
  archive & ifChooseNug;
  archive & ifAssumeRcondZero;
  archive & ifPrescribedNug;
  archive & nug;
  archive & iPtsKeep;
  archive & numPoints;
  archive & numPointsKeep;
  archive & numWholePointsKeep;
  archive & numExtraDerKeep;
  archive & numEqnAvail;
  archive & numRowsR;
  archive & ifHaveAnchorPoint;
  archive & iAnchorPoint;
  //don't archive XR since it's only a reference into SurfData .xr
  archive & XRreorder;
  //archive & Yall; //need this during the construction of a model but not afterward so don't archive
  archive & Y;
  //don't archive Gall we need this during the construction of a model but not afterward
  //don't archive Gtran we need this during the construction of a model but not afterward
  archive & ifReducedPoly;
  archive & polyOrderRequested;
  archive & polyOrder;
  archive & numTrend;
  archive & nTrend;
  //don't archive iTrendKeep, it's not needed because at the discarded terms in the trend basis function are removed from Poly at the end of create()
  archive & Poly;
  //don't archive flyPoly, it's work space needed to efficiently evaluate a polynomial basis "on the fly" (it's a member variable only to avoid constant allocation and deallocation)
  archive & betaHat;
  //don't archive derivBeta, it's work space needed to efficiently evaluate derivatives of a polynomial (it's a member variable only to avoid constant allocation and deallocation)
  //don't archive Z, we need it during the construction of a model but not afterward
  //don't archive Ztran_theta, we need it during the construction of a model but not afterward
  //don't archive deltaXR, we need it during the construction of a model but not afterward
  //don't archive R, we need it during the construction of a model but not afterward, once we have ranking of candidate points by pivoteted cholesky we will use it as temporary variable space after the model is constructed but it still won't be something we want to retain
  //archive & Rinv; //not used, would be used for analytic integral of adjusted variance
  archive & RChol;
  //don't archive scaleRChol, we need it during the construction of a model but not afterward
  //don't archive sumAbsColR, we need it during the construction of a model but not afterward
  //don't archive oneNormR, we need it during the construction of a model but not afterward
  //don't archive lapackRcondR, weneed it during the construction of a model but not afterward
  //don't archive rcondDblWork, we need it during the construction of a model but not afterward
  //don't archive rcondIntWork, we need it during the construction of a model but not afterward
  archive & rcondR;
  archive & rcond_G_Rinv_Gtran;
  archive & Rinv_Gtran;
  //don't archive G_Rinv_Gtran, we need it during the construction of a model but not afterward
  archive & G_Rinv_Gtran_Chol;
  //don't archive G_Rinv_Gtran_Chol_Scale, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Gtran_Chol_DblWork, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Gtran_Chol_IntWork, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Y, we need it during the construction of a model but not afterward
  //don't archive eps, we need it during the construction of a model but not afterward
  archive & rhs;
  archive & estVarianceMLE;
  archive & likelihood;
  //don't archive prevObjDerMode, we need it during the construction of a model but not afterward
  //don't archive prevConDerMode, we need it during the construction of a model but not afterward
  //don't archive prevTheta, we need it during the construction of a model but not afterward
  archive & maxObjDerMode;
  archive & maxConDerMode;
  archive & obj;
  //don't archive con, we need it during the construction of a model but not afterward
}
#endif

#endif
