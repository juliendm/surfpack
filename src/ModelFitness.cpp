/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "SurfData.h"
#include "SurfpackModel.h"
#include "ModelFitness.h"
#include "ModelFactory.h"

using std::cout;
using std::endl;
using std::max_element;
using std::min_element;
using std::random_shuffle;
using std::set;
using std::string;
using std::vector;
using surfpack::shared_rng;


// --------------------------
// implementation of Residual
// --------------------------

Residual::Residual(DifferenceType dt_in) : dt(dt_in) 
{ /* empty ctor */ }


double Residual::operator()(double observed, double predicted) const
{
    switch(dt) {
      case DT_ABSOLUTE: return fabs(observed - predicted);
      case DT_SCALED: return fabs(observed - predicted)/fabs(observed);
      case DT_SQUARED: return (observed-predicted)*(observed-predicted);
      default: assert(false);
    }
    assert(dt == DT_ABSOLUTE || dt == DT_SCALED || dt == DT_SQUARED); 
    return 0.0;
}


// ----------------------------
// implementation of VecSummary
// ----------------------------

VecSummary::VecSummary(MetricType mt_in) 
  : mt(mt_in) 
{ /* empty ctor */ }


double VecSummary::operator()(const VecDbl& resids) const
{
  switch (mt) {
    case MT_SUM: return std::accumulate(resids.begin(),resids.end(),0.0);
    case MT_MEAN: return surfpack::mean(resids);
    case MT_ROOT_MEAN: return sqrt(surfpack::mean(resids));
    case MT_MAXIMUM: 
      VecDbl::const_iterator itr = max_element(resids.begin(),resids.end());
      return *itr;
    //default: throw string("Unknown vec summary");
  }
  return 0.0;
}


// ------------------------------
// implementation of ModelFitness
// ------------------------------

double ModelFitness::operator()(const VecDbl& obs, const VecDbl& pred) const
{
  throw string("Not implemented for abstract ModelFitness class");
}

ModelFitness* ModelFitness::Create(const std::string& metric, unsigned n)
{
  if (metric == "sum_squared") {
    return new StandardFitness(Residual(DT_SQUARED),VecSummary(MT_SUM));
  } else if (metric == "mean_squared") {
    return new StandardFitness(Residual(DT_SQUARED),VecSummary(MT_MEAN));
  } else if (metric == "root_mean_squared") {
    return new StandardFitness(Residual(DT_SQUARED),VecSummary(MT_ROOT_MEAN));
  } else if (metric == "max_squared") {
    return new StandardFitness(Residual(DT_SQUARED),VecSummary(MT_MAXIMUM));
  } else if (metric == "sum_scaled") {
    return new StandardFitness(Residual(DT_SCALED),VecSummary(MT_SUM));
  } else if (metric == "mean_scaled") {
    return new StandardFitness(Residual(DT_SCALED),VecSummary(MT_MEAN));
  } else if (metric == "max_scaled") {
    return new StandardFitness(Residual(DT_SCALED),VecSummary(MT_MAXIMUM));
  } else if (metric == "sum_abs") {
    return new StandardFitness(Residual(DT_ABSOLUTE),VecSummary(MT_SUM));
  } else if (metric == "mean_abs") {
    return new StandardFitness(Residual(DT_ABSOLUTE),VecSummary(MT_MEAN));
  } else if (metric == "max_abs") {
    return new StandardFitness(Residual(DT_ABSOLUTE),VecSummary(MT_MAXIMUM));
  } else if (metric == "press") {
    return new PRESSFitness();
  } else if (metric == "cv") {
    return new CrossValidationFitness(n);
  } else if (metric == "rsquared") {
    return new R2Fitness();
  }
  string msg = "Metric '" + metric + "' not supported";
  throw msg; 
  return new StandardFitness(Residual(DT_SQUARED),VecSummary(MT_SUM));
}


VecDbl ModelFitness::getResiduals(const Residual& resid, 
				  const VecDbl& obs, const VecDbl& pred)
{
  assert(obs.size() == pred.size());
  VecDbl result(obs.size());
  for (unsigned i = 0; i < result.size(); i++) {
    result[i] = resid(obs[i],pred[i]);
  }
  return result;
}


// ---------------------------------
// implementation of StandardFitness
// ---------------------------------

StandardFitness::StandardFitness()
: resid(Residual(DT_SQUARED)), vecsumry(MT_MEAN)
{ /* empty ctor */ }


StandardFitness::StandardFitness(const Residual& resid_in, 
				 const VecSummary& vecsumry_in)
: resid(resid_in), vecsumry(vecsumry_in)
{ /* empty ctor */ }


double StandardFitness::operator()(const VecDbl& obs, const VecDbl& pred) const
{
  VecDbl residuals = getResiduals(resid,obs,pred);
  return vecsumry(residuals);
}


double StandardFitness::operator()(const SurfpackModel& sm, 
				   const SurfData& sd) const
{
  VecDbl predicted = sm(sd);
  VecDbl observed = sd.getResponses();
  VecDbl residuals = getResiduals(resid,observed,predicted);
  return vecsumry(residuals);
}


// ----------------------------------------
// implementation of CrossValidationFitness
// ----------------------------------------

// BMA TODO: consider moving to root_mean_squared for default metric
CrossValidationFitness::CrossValidationFitness()
  : ModelFitness(), num_partitions(10), default_metric("mean_squared")
{ /* empty ctor */ }


// BMA TODO: consider moving to root_mean_squared for default metric
CrossValidationFitness::CrossValidationFitness(unsigned n_in)
  : ModelFitness(), num_partitions(n_in), default_metric("mean_squared")
{ /* empty ctor */ }


/** Historically CV used a single mean-squared fitness metric */
double CrossValidationFitness::
operator()(const SurfpackModel& sm, const SurfData& sd) const
{
  // get predictions (estimates) from the cross-validation core
  VecDbl estimates;
  leaveout_estimates(estimates, sm, sd);
  //cout << "CV vals: " << surfpack::fromVec<double>(estimates) << endl;

  // get the active response from SurfData
  VecDbl responses = sd.getResponses();

  return calc_one_metric(responses, estimates, default_metric);
}


/** Perform a single set of CV leave-out builds, but compute multiple metrics */
void CrossValidationFitness::
eval_metrics(VecDbl& metric_values, const SurfpackModel& sm, const SurfData& sd,
	     const VecStr& metric_names) const
{
  // get predictions (estimates) from the cross-validation core
  VecDbl estimates;
  leaveout_estimates(estimates, sm, sd);

  // get the active response from SurfData
  VecDbl responses = sd.getResponses();

  // now compute all the metrics requested
  metric_values.clear();
  metric_values.reserve(metric_names.size());
  for (VecStr::const_iterator mi = metric_names.begin();
       mi != metric_names.end(); ++ mi)
    metric_values.push_back(calc_one_metric(responses, estimates, *mi));
}


void CrossValidationFitness::
leaveout_estimates(VecDbl& estimates, const SurfpackModel& sm,
		   const SurfData& sd) const
{
  // The data is divided into n partitions for cross validation.  
  unsigned n_final = num_partitions;

  // When n = 0, leave out about 10% of the data (n=10 partitions),
  // but with a upper bound of the data size (results in PRESS).
  const unsigned default_partitions = 10;
  if (num_partitions == 0)
    n_final = default_partitions;

  // enforce a lower bound of 2 and upper bound of data size
  const unsigned min_partitions = 2;
  if (n_final > sd.size())
    n_final = sd.size();
  else if (n_final < min_partitions)
    n_final = min_partitions;

  /*if you want cross validation leaving m out then n=sd.size()/m

    low  = partition*my_data.size()/n
    high = (partition+1)*my_data.size()/n-1

    partition loops from 0 to n-1 

    say partition=0   => low =0*(n+1)/n=0
                         high=1*(n+1)/n-1=1+1/n-1 = 0
    say partition=n-1 => low =(n-1)*(n+1)/n = (n^2-1)/n=n
                         high =n*(n+1)/n-1= n+1-1 = n

    if n=my_data.size()-2 then
    say partition=0   => low =0*(n+2)/n  = 0
                         high=1*(n+2)/n-1= 0

    if n=my_data.size()/2 then
    say partition=0   => low =0*2*n/n    = 0
                         high=1*2*n/n-1  = 1
    say partition=n-1 => low =(n-1)*2*n/n= 2*n-2
                         high= n*2*n/n-1 = 2*n-1

    if n=my_data.size()/3 then
    say partition=0   => low =0
                      => high=1*3*n/n-1=2
    say partition=n-1 => low =(n-1)*3*n/n= 3*n-3
                      => high=n*3*n/n-1  = 3*n-1
  */

  //cout << "CV Fitness: " << n << endl;
  SurfData my_data = sd; // Get non const copy
  ParamMap args = sm.parameters();

  // silence model output for cross-validation
  args["verbosity"] = surfpack::toString<short>(surfpack::SILENT_OUTPUT);

  VecUns indices(my_data.size()); 
  for (unsigned i = 0; i < indices.size(); i++) indices[i] = i;
  random_shuffle(indices.begin(),indices.end(),shared_rng());
  estimates.resize(my_data.size());
  for (unsigned partition = 0; partition < n_final; partition++) {
    //cout << "part: " << partition << endl;
    SetUns excludedPoints;
    unsigned low = surfpack::block_low(partition, n_final, my_data.size());
    unsigned high = surfpack::block_high(partition, n_final, my_data.size());
    //cout << "low/high: " << low << " " << high << endl;
    for (unsigned k = low; k <= high; k++) excludedPoints.insert(indices[k]);
    my_data.setExcludedPoints(excludedPoints);
    //cout << " excludes: " << excludedPoints.size() << endl;
    SurfpackModelFactory* factory = ModelFactory::createModelFactory(args);
    SurfpackModel* model = factory->Build(my_data);
    my_data.setExcludedPoints(SetUns());
    for (unsigned k = low; k <= high; k++) {
      estimates[indices[k]] = (*model)(my_data(indices[k]));
      //cout << "for k = " << k << ": " << estimates[indices[k]] << endl;
    }
    delete model;
    delete factory;
  }
}


double CrossValidationFitness::
calc_one_metric(const VecDbl& observed, const VecDbl& predicted,
		const string& metric_name) const
{
  assert(observed.size() == predicted.size());

  ModelFitness* mf = ModelFitness::Create(metric_name);
  double fitness = (*mf)(predicted, observed);
  delete mf;

  return fitness;
}


// ------------------------------
// implementation of PRESSFitness
// ------------------------------

PRESSFitness::PRESSFitness()
{ /* empty ctor */ }


double PRESSFitness::operator()(const SurfpackModel& sm, 
				const SurfData& sd) const
{
  // create a leave one out cross-validation fitness operator here
  // since the data size isn't available at construct time (number of
  // partitions equals number of points) TODO: for efficiency, might
  // want to just reimplement instead of calling CV fitness; for now,
  // assume model construction dominates bookeeping arithmetic
  // overhead
  ModelFitness* cvmf = ModelFitness::Create("cv", sd.size());
  double fitness = (*cvmf)(sm, sd);
  delete cvmf;
  return fitness;
}


// ---------------------------------
// implementation of R2Fitness
// ---------------------------------

R2Fitness::R2Fitness()
{ /* empty ctor */ }


double R2Fitness::operator()(const SurfpackModel& sm, const SurfData& sd) const
{

  VecDbl predicted = sm(sd);
  VecDbl observed = sd.getResponses();

  return this->operator()(observed, predicted);
}


double R2Fitness::operator()(const VecDbl& obs, const VecDbl& pred) const
{
  double obs_mean = surfpack::mean(obs);
  VecDbl vec_mean = VecDbl(obs.size(),obs_mean);
  StandardFitness sum_squares = 
    StandardFitness(Residual(DT_SQUARED),VecSummary(MT_SUM));

  return sum_squares(pred,vec_mean)/sum_squares(obs,vec_mean);
}
