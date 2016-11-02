/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_MODEL_H__
#define __SURFPACK_MODEL_H__

#include "surfpack_system_headers.h"
#include "surfpack.h"
// need complete type information for serialization
#include "ModelScaler.h"

class SurfData;


///////////////////////////////////////////////////////////
///	Surfpack Model 
///////////////////////////////////////////////////////////

class SurfpackModel
{

public:

  SurfpackModel(unsigned ndims_in);
  SurfpackModel(const SurfpackModel& other);
  virtual VecDbl operator()(const SurfData& data) const;
  double operator()(const VecDbl& x) const;
  virtual double variance(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual MtxDbl hessian(const VecDbl& x) const;
  virtual std::string asString() const = 0;
  /// return the value of some error metric
  double goodnessOfFit(const std::string MetricName, const SurfData& data); 
  ///double goodnessOfFit(const std::string MetricName); 
  /// Return R^2, which measures the proportion of variability in the data 
  /// accounted for by the model (the approximating surface).
  /// R^2 = 1-SSE/SST = SSR/SST. 
  double rSquared(const SurfData& data);
  /// For each point x in dataSet, construct a similar approximation Surface
  /// that includes all of the points in dataSet except x.  Then evaluate the
  /// Surface at x.  The difference between x and the estimate of x given the
  /// rest of the data is the residual for x.  The PRESS statistic is the
  /// square root of the mean of the squares of all the residuals.
  double press(const SurfData& data);
  double nFoldCrossValidation(const SurfData& data, unsigned n);
  /// Compute one of several goodness of fit metrics.  The observed parameter
  /// should be a list of observed (or true) function values; the vector of
  /// predicted values gives the corresponding estimates from this surface.
  /// The dt parameter specifies the kind of residuals to compute.  DT_ABSOLUTE
  /// residuals are (observed - predicted), DT_SQUARED residuals are the squares
  /// of the absolute residuals.  DT_SCALED residuals are the ABSOLUTE residuals
  /// divided by the observed value.  Given the type of residuals, the client
  /// may request the min, max, sum, or mean of the set of residuals over all
  /// the given data points.  Two additional metrics are possible.  The
  /// relative maximum absolute error is the maximum absolute error divided
  /// by the standard deviation of the observed data.  The relative average
  /// absolute error is the mean absolute error divided by the standard
  /// deviation of observed data.
  double genericMetric(std::vector<double>& observed,
    std::vector<double>& predicted, enum MetricType mt, enum DifferenceType dt);
  virtual ~SurfpackModel();
  virtual void scaler(ModelScaler* ms);
  ModelScaler* scaler() const;
  unsigned size() const { return ndims;}
  const ParamMap& parameters() const { return args; }
  void parameters(const ParamMap& args) { this->args = args;}

protected:

  /// default constructor used when reading from archive file 
  SurfpackModel();

  /// evaluation function implemented by derived classes, used in operator()
  virtual double evaluate(const VecDbl& x) const = 0;

  /// number of input (x) variables
  unsigned ndims;
  /// model configuration parameters
  ParamMap args;
  /// data scaler for this model
  ModelScaler*  mScaler;

private:

  /// disallow assignment as not implemented
  SurfpackModel& operator=(const SurfpackModel& other);

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for base class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


///////////////////////////////////////////////////////////
///	Surfpack Model Factory
///////////////////////////////////////////////////////////

class SurfpackModelFactory
{

public:

  /// Default constructor 
  SurfpackModelFactory();
  /// Partial parameter list constructor
  SurfpackModelFactory(const ParamMap& args);

  /// Build a model from the provided SurfData
  virtual SurfpackModel* Build(const SurfData& sd);

  /// the minimum number of points with which Surfpack will build a model
  virtual unsigned minPointsRequired();
  /// the recommended default number of points
  virtual unsigned recommendedNumPoints();
  /// whether the model type supports constraints (anchor point)
  virtual bool supports_constraints();

  /// retreive the configuration parameters
  const ParamMap& parameters() const;
  /// add a configuration parameter
  void add(const std::string& name, const std::string& value);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd) = 0;

  /// set member data prior to build; derived classes are responsible
  /// for calling this implementation
  virtual void config();

  /// convenience function to verify that a model has sufficient data to build
  virtual void sufficient_data(const SurfData& sd);

  /// map of configuration parameters
  ParamMap params;
  /// dimension of the problem (variables)
  unsigned ndims;
  /// active response index over which to build
  unsigned response_index;

};


// -----
// Definitions and export of serialization functions
// -----

/** Serializer for the base class data, e.g., parameter list.  This is
    called by the derived class serialize functions via base_object */
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template<class Archive> 
void SurfpackModel::serialize(Archive & archive, const unsigned int version)
{
  archive & args;
  archive & ndims;
  archive & mScaler;
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(SurfpackModel)

#endif


#endif  // __SURFPACK_MODEL_H__
