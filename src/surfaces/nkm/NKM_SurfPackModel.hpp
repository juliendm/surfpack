#ifndef _SURFPACKMODEL_H_
#define _SURFPACKMODEL_H_

#include "NKM_SurfMat.hpp"
#include "NKM_SurfPack.hpp"
#include "NKM_Optimize.hpp"
#include <cfloat>
#include <iostream>
#include <exception>

namespace nkm {

class SurfPackModel 
{
protected:
  SurfData sdBuild;  
  SurfDataScaler scaler;
  short outputLevel;

private:
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for base class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

public:
  
  SurfPackModel() : scaler(sdBuild), outputLevel(NORMAL_OUTPUT) {};

  SurfPackModel(const SurfData& sd,int iout_keep) : sdBuild(sd,iout_keep), scaler(sdBuild), outputLevel(NORMAL_OUTPUT) {};

  virtual void create() {
    std::cerr << "the create() function has not been implemented for this model type" << std::endl;
    return;
  };

  virtual std::string model_summary_string() const {
    std::string mod_sum_str="the model_summary_string() function has not been implemented for this model\n";
    return mod_sum_str;
  };


  virtual double evaluate(const  MtxDbl& xr) =0;

  virtual MtxDbl& evaluate(MtxDbl& y, const MtxDbl& xr) 
  {
    int nvarsxr=xr.getNRows();
    int nptsxr=xr.getNCols();
    assert((nvarsxr==sdBuild.getNVarsr())&&(nptsxr>0));
    y.newSize(1,nptsxr);

    if(nptsxr==1) {
      y(0,0)=evaluate(xr);
      return y;
    }
      
    MtxDbl xr_temp(nvarsxr,1);
    for(int ipt=0; ipt<nptsxr; ++ipt) {
      xr.getCols(xr_temp,ipt);
      y(0,ipt)=evaluate(xr_temp);
    }
    return y;
  };

  virtual double eval_variance(const MtxDbl& xr) {
    std::cerr << "This model doesn't have an implemented function to return a variance" << std::endl;
    assert(false);
    // stricter compilers don't allow divide by 0
    // consider use of cmath macro NAN, though might not port to MSVS
    //return (0.0/0.0);
    // since code unreachable anyway, using a large variance value
    return(DBL_MAX);
  };

  virtual MtxDbl& eval_variance(MtxDbl& var, const MtxDbl& xr) 
  {
    int nvarsxr=xr.getNRows();
    int nptsxr=xr.getNCols();
    assert((nvarsxr==sdBuild.getNVarsr())&&(nptsxr>0));
    var.newSize(1,nptsxr);

    if(nptsxr==1) {
      var(0,0)=eval_variance(xr);
      return var;
    }
      
    MtxDbl xr_temp(nvarsxr,1);
    for(int ipt=0; ipt<nptsxr; ++ipt) {
      xr.getCols(xr_temp,ipt);
      var(0,ipt)=eval_variance(xr_temp);
    }
    return var;
  };

  virtual MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) =0;

  virtual MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) =0;


  /// adjust correlations to be feasible with respect to condition
  /// number constraints
  virtual MtxDbl& makeGuessFeasible(MtxDbl& correlations, 
				    OptimizationProblem *opt)
  {
    // default at base class is no-op
    return correlations;
  };

  virtual void getRandGuess(MtxDbl& guess) const{


  };

  virtual void set_conmin_parameters(OptimizationProblem& opt) const{
  };

  virtual void set_direct_parameters(OptimizationProblem& opt) const{
  };
  

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "good" KrigingModel)
  virtual double objective(const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objective" << std::endl;
    throw(std::string("Derived does not implement"));
  }

  /// the objective function, i.e. the negative log(likelihood), and
  /// its gradient; minimizing the objective function produces a good
  /// KrigingModel
  virtual void objectiveAndGradient(double& Obj, MtxDbl& GradObj,
				    const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndGradient" << std::endl;
    throw(std::string("Derived does not implement"));
  }  


  /// objective plus condition number constraints
  virtual void objectiveAndConstraints(double& Obj, MtxDbl& Con, 
				       const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndConstraints" << std::endl;
    throw(std::string("Derived does not implement"));
  }


  /// objective plus condition number constraints with gradients
  virtual void objectiveAndConstraintsAndGradients(double& Obj, MtxDbl& Con, 
						   MtxDbl& GradObj, 
						   MtxDbl& GradCon, 
						   const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndConstraintsAndGradients" << std::endl;
    throw(std::string("Derived does not implement"));
  }



};

} // end namespace nkm

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template<class Archive> 
void nkm::SurfPackModel::serialize(Archive & archive, 
				   const unsigned int version)
{
  archive & sdBuild;
  archive & scaler;
  archive & outputLevel;
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(nkm::SurfPackModel)
#endif

#endif
