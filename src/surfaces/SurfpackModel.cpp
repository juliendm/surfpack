#include "SurfpackModel.h"
#include "SurfData.h"
#include "surfpack.h"
#include "ModelScaler.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(SurfpackModel)
#endif

///////////////////////////////////////////////////////////
///	Surfpack Model 
///////////////////////////////////////////////////////////

VecDbl SurfpackModel::operator()(const SurfData& data) const
{
  VecDbl result(data.size());
  for (unsigned pt = 0; pt < data.size(); pt++) {
    result[pt] = (*this)(data(pt));
  }
  return result;
}

double SurfpackModel::operator()(const VecDbl& x) const
{
  //cout << "\nunscaled\n";
  //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout," "));
  const VecDbl& x1 = mScaler->scale(x);
  //cout << "\nscaled\n";
  //copy(x1.begin(),x1.end(),std::ostream_iterator<double>(cout," "));
  double value = evaluate(x1);
  //cout << "evaluated: " << value << "\n";
  double result = mScaler->descale(value);
  //cout << "descaled: " << result << endl;
  return result;
  
  //return mScaler->descale(evaluate(mScaler->scale(x)));
}


SurfpackModel::SurfpackModel(): 
  ndims(0), mScaler(NULL)
{ /* empty ctor */ }

SurfpackModel::SurfpackModel(unsigned ndims_in) 
  : ndims(ndims_in), mScaler(new NonScaler)
{

}

SurfpackModel::SurfpackModel(const SurfpackModel& other)
  : ndims(other.ndims), args(other.args), mScaler(other.mScaler->clone())
{

}

SurfpackModel::~SurfpackModel()
{
  if (mScaler) {
    delete mScaler;
    mScaler = NULL;
  }
}

double SurfpackModel::variance(const VecDbl& x) const
{
  throw std::string("This model does not currently support variance eval");
}

VecDbl SurfpackModel::gradient(const VecDbl& x) const
{
  throw std::string("This model does not currently support gradients");
}

MtxDbl SurfpackModel::hessian(const VecDbl& x) const
{
  throw std::string("This model does not currently support hessians");
}

double SurfpackModel::goodnessOfFit(const string metricName, const SurfData& surf_data) 
{
  // need to put error message if surfData not populated
  std::cout << " Calculating goodness of fit " << "\n" ;
  if (metricName == "rSquared") {
    return rSquared(surf_data);
  } else if (metricName == "press") {
    return nFoldCrossValidation(surf_data,surf_data.size());
  } else {
    // The rest of these metrics all have many computations in common
    // and are grouped together in the genericMetric method
    vector<double> observed(surf_data.size());
    vector<double> predicted(surf_data.size());
    observed = surf_data.getResponses();
    for (unsigned i = 0; i < surf_data.size(); i++) {
     predicted[i]=(*this)(surf_data(i));
    }
    //predicted = surf_data.getPredictor();
    //getValue(surfData, observed, predicted);
    if (metricName == "min_abs" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,DT_ABSOLUTE);
    } else if (metricName == "max_abs") {
      return genericMetric(observed,predicted,MT_MAXIMUM,DT_ABSOLUTE);
    } else if (metricName == "sum_abs") {
      return genericMetric(observed,predicted,MT_SUM,DT_ABSOLUTE);
    } else if (metricName == "mean_abs") {
      return genericMetric(observed,predicted,MT_MEAN,DT_ABSOLUTE);
    } else if (metricName == "max_relative") {
      return genericMetric(observed,predicted,MT_RELATIVE_MAXIMUM,DT_ABSOLUTE);
    } else if (metricName == "mean_relative") {
      return genericMetric(observed,predicted,MT_RELATIVE_AVERAGE,DT_ABSOLUTE);
    } else if (metricName == "min_squared" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,DT_SQUARED);
    } else if (metricName == "max_squared") {
      return genericMetric(observed,predicted,MT_MAXIMUM,DT_SQUARED);
    } else if (metricName == "sum_squared") {
      return genericMetric(observed,predicted,MT_SUM,DT_SQUARED);
    } else if (metricName == "mean_squared") {
      return genericMetric(observed,predicted,MT_MEAN,DT_SQUARED);
    } else if (metricName == "min_scaled" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,DT_SCALED);
    } else if (metricName == "max_scaled") {
      return genericMetric(observed,predicted,MT_MAXIMUM,DT_SCALED);
    } else if (metricName == "sum_scaled") {
      return genericMetric(observed,predicted,MT_SUM,DT_SCALED);
    } else if (metricName == "mean_scaled") {
      return genericMetric(observed,predicted,MT_MEAN,DT_SCALED);
    //} else if (metricName == "root_mean_squared") {
    //  return rootMeanSquared(observed,predicted);
    //} else if (metricName == "l2norm") {
    //  return sqrt(genericMetric(observed,predicted,MT_SUM,DT_SQUARED));
    } else {
      throw string("No error metric of that type in this class");
    }
  }
}
 
double SurfpackModel::rSquared(const SurfData& surf_data)
{
  // Sum of the function evaluations for all the points; used to compute mean
  double sumObserved = 0.0;
 
  // Sum of the squared function evaluations for the points; used to compute
  // the total sum of squares
  double sumOfSquaresObserved = 0.0;
 
  // Sum of the squared differences between the true function value and the
  // Surface approximation's estimate of the function value over all of the
  // points
  double residualSumOfSquares = 0.0;
 
  // Sum of the squared differences between the true function value and mean
  // function value over all of the points
  double totalSumOfSquares = 0.0;
 
  for (unsigned i = 0; i < surf_data.size(); i++) {
    double observedF = surf_data.getResponse(i);
    std::cout << "observedF " << observedF;
    double estimatedF = (*this)(surf_data(i));
    std::cout << "estimatedF " << estimatedF;
    //double estimatedF = surf_data.getPredictor(i);
    double residual = observedF - estimatedF;
    std::cout << "residual " << residual << "\n";
    residualSumOfSquares += residual * residual;
    sumObserved += observedF;
    sumOfSquaresObserved += observedF * observedF;
  }
  // This is the same as sigma{i=1 to n}(x_i - xbar)^2
  totalSumOfSquares = sumOfSquaresObserved -
    (sumObserved * sumObserved / surf_data.size());
  double rSquaredValue = 1.0 - residualSumOfSquares / totalSumOfSquares;
  // In a polynomial regression, this computation will always result in a value
  // between 0 and 1, because the residual sum of squares will always be less
  // than or equal to the total sum of squares (i.e., the worst your regression
  // can do is give you the mean everywhere).  Some of the other methods (Mars,
  // Kriging, RBF, etc.) can have much larger residual sum of squares.  We will
  // refrain from returning a value less than zero however.
  return (rSquaredValue < 0) ? 0 : rSquaredValue;
}

double SurfpackModel::press(const SurfData& data)
{
return 0.;
}

double SurfpackModel::nFoldCrossValidation(const SurfData& data, unsigned n)
{
return 0.;
}

double SurfpackModel::genericMetric(std::vector<double>& observed,
    std::vector<double>& predicted, enum MetricType mt, enum DifferenceType dt)
{
return 0.;
}

void SurfpackModel::scaler(ModelScaler* ms)
{
  delete mScaler;
  mScaler = ms->clone();
}

ModelScaler* SurfpackModel::scaler() const
{
  return mScaler;
}


///////////////////////////////////////////////////////////
///	Surfpack Model Factory
///////////////////////////////////////////////////////////

/** Default constructor for a factory that will get its configuration
    later, likely from a SurfData object */
SurfpackModelFactory::SurfpackModelFactory()
  : params(), ndims(0), response_index(0)
{

}

/** Factory constructor that has some of its parameters passed
    initially, but not necessarily the number of dimensions */
SurfpackModelFactory::SurfpackModelFactory(const ParamMap& params_in)
  : params(params_in), ndims(0), response_index(0)
{

}

/** Use information in the parameter map to set member data prior to
    build/query.  Derived classes are responsible for calling this
    implementation. */
void SurfpackModelFactory::config()
{
  ndims = std::atoi(params["ndims"].c_str());
  assert(ndims);
  string arg = params["response_index"];
  if (arg != "") response_index = std::atoi(arg.c_str());
}

/// Default implementation of minimum data to build
unsigned SurfpackModelFactory::minPointsRequired()
{
  SurfpackModelFactory::config();
  assert(ndims);
  return (ndims+1);
}

/// Default implementation of recommended data to build
unsigned SurfpackModelFactory::recommendedNumPoints()
{
  SurfpackModelFactory::config();
  assert(ndims);
  return (5*ndims);
}

/// By default, models don't support constraints
bool SurfpackModelFactory::supports_constraints()
{
  return false;
}


// default implementation assumes only function value data
void SurfpackModelFactory::sufficient_data(const SurfData& sd)
{
  if (sd.size() < minPointsRequired()) {
    std::ostringstream not_enough;
    not_enough << "Not enough Points: ";
    not_enough << "size of data = " << sd.size();
    not_enough << ", minPointsRequired = " << minPointsRequired();
    throw(not_enough.str());
  }
}

const ParamMap& SurfpackModelFactory::parameters() const
{
  return params;
}

void SurfpackModelFactory::add(const std::string& name, const std::string& value)
{
  params[name]=value;
}

/** Update the model factory with characteristics from the provided
    data, select the desired response index in the data, and create a
    SurfpackModel. Hence the same factory can be used to build models
    on different responses */
SurfpackModel* SurfpackModelFactory::Build(const SurfData& sd)
{
  //cout << "Data:\n";
  //sd.writeText(cout);
  //cout << "\nParams:\n";
  //
  //  for (ParamMap::iterator itr = params.begin();
  //      itr != params.end(); itr++) {
  //    std::cout << "     " << itr->first << ": " << itr->second << std::endl;
  //  }
  this->add("ndims",surfpack::toString<unsigned>(sd.xSize()));
  this->config();
  // TODO: add verbosity control
  if (false) {
    std::cout << "SurfpackModelFactory built with parameters:" << std::endl; 
    for (ParamMap::iterator itr = params.begin();
	 itr != params.end(); itr++) {
      std::cout << "     " << itr->first << ": " << itr->second << std::endl;
    }
  }
  sd.setDefaultIndex(this->response_index);
  // check whether there is sufficient SurfData to build
  sufficient_data(sd);
  SurfpackModel* model = Create(sd);
  model->parameters(params);
  return model;
}
