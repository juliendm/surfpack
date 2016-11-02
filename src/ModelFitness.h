/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MODEL_FITNESS_H__
#define __MODEL_FITNESS_H__

#include "surfpack_system_headers.h"
#include "surfpack.h"

/// class to compute a single scalar residual, with a specific difference type
class Residual
{
public:

  /// standard constructor for Residual, accepting a DifferenceType
  Residual(DifferenceType dt_in);

  /// evaluate the residual given observed and predicted
  double operator()(double observed, double predicted) const; 

protected:

  /// type of residual difference as defined in surfpack.h
  DifferenceType dt;

private:

  /// make the default constructor private, since there's no default diff type
  Residual();
};


/// class VecSummary to compute a summary statistic over a vector of
/// residuals; used in StandardFitness
class VecSummary
{
public:

  /// standard constructor for VecSummary, accepting a MetricType
  VecSummary(MetricType mt_in);

  /// compute the configured vector summary for the passed redisuals 
  double operator()(const VecDbl& resids) const;

protected:

  /// type of vector metric as defined in surfpack.h
  MetricType mt;

private:

  /// make the default constructor private, since there's no default diff type
  VecSummary();
};

/// Generic ModelFitness assessor
class ModelFitness 
{
public:

  virtual double operator()(const SurfpackModel& sm, 
			    const SurfData& sd) const = 0;

  // BMA: Do we need this? Would be better to have a pure virtual base
  
  virtual double operator()(const VecDbl& obs, const VecDbl& pred) const;

  /// factory to return a derived ModelFitness type by name
  /// n is used to configure CrossValidationFitness folds
  static ModelFitness* Create(const std::string& metric, unsigned n = 0);

  /// get residuals, transformed by the passed Residual modifier resid
  static VecDbl getResiduals(const Residual& resid, 
			     const VecDbl& obs, const VecDbl& pred);
};


/// StandardFitness applies a Residual type to each residual and then
/// computes the VecSummary of them
class StandardFitness : public ModelFitness
{
public:

  /// default constructor using mean squared residuals
  StandardFitness();

  /// standard constructor accepting a residual and summary tupe
  StandardFitness(const Residual& resid_in, const VecSummary& vecsumry_in);

  /// compute the fitness for the passed model and active response in sd
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;

  /// compute the fitness for the passed list of observations and predictions
  virtual double operator()(const VecDbl& obs, const VecDbl& pred) const;

protected:

  /// residual operator to apply to each difference
  Residual resid;

  /// summary operator to apply over all data points
  VecSummary vecsumry;
};


/// k-fold cross validation fitness
/** Partition data into num_partitions partitions.  For each, rebuild
    the model leaving out the partition, compute residuals against the
    leave out data. */
class CrossValidationFitness : public ModelFitness
{
public:
  
  /// default constructor, setting up default folds = 10
  CrossValidationFitness();

  /// constructor accepting number of cross-validation partitions
  CrossValidationFitness(unsigned n_in);

  /// compute the fitness for the passed model and active response in sd
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;

  /// perform cross-validation, returning statistics for multiple metrics
  void eval_metrics(VecDbl& metrics, const SurfpackModel& sm,
		    const SurfData& sd, const VecStr& metric_names) const;

protected:

  /// core to calculate predictions leaving out a sequence of folds (partitions)
  void leaveout_estimates(VecDbl& estimates, const SurfpackModel& sm,
			  const SurfData& sd) const;

  /// calculate a single fitness metric for the cross validation data
  double calc_one_metric(const VecDbl& observed, const VecDbl& predicted,
			 const std::string& metric_name) const;

  /// number of partitions (folds) for cross-validation
  unsigned num_partitions;

  /// default metric to use in operator(); historically mean_squared
  std::string default_metric;
  
};


/// leave one out cross validation fitness; perhaps a misnomer to call PRESS
class PRESSFitness: public ModelFitness
{
public:
  PRESSFitness();
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
};


/// polynomial regression R-squared fitness: prediction variance
/// w.r.t. observation mean divided by observed variance
/// w.r.t. observation mean
class R2Fitness: public ModelFitness
{
public:
  R2Fitness();
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
  virtual double operator()(const VecDbl& obs, const VecDbl& pred) const;
};

#endif
