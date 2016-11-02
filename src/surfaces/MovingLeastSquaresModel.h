/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MOVING_LEAST_SQUARES_MODEL_H__
#define __MOVING_LEAST_SQUARES_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfData.h"
#include "LinearRegressionModel.h"

class MovingLeastSquaresModel : public SurfpackModel
{
public:
  MovingLeastSquaresModel(const SurfData& sd_in, const LRMBasisSet& bs_in,
    unsigned continuity_in = 1);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  SurfData sd;
  LRMBasisSet bs;
  mutable VecDbl coeffs;
  unsigned continuity;
  
friend class MovingLeastSquaresModelTest;

private:
  /// default constructor used when reading from archive file 
  MovingLeastSquaresModel() { /* empty ctor */}

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

  /// disallow copy construction as not implemented
  MovingLeastSquaresModel(const MovingLeastSquaresModel& other);

  /// disallow assignment as not implemented
  MovingLeastSquaresModel& operator=(const MovingLeastSquaresModel& other);

};

///////////////////////////////////////////////////////////
///	Moving Least Squares Model Factory
///////////////////////////////////////////////////////////

class MovingLeastSquaresModelFactory : public SurfpackModelFactory 
{

public:
  MovingLeastSquaresModelFactory();
  MovingLeastSquaresModelFactory(const ParamMap& args);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  unsigned weight;
  unsigned order;
};

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
/** Serializer for the dervied Model data, e.g., basis and coefficients.
    Must call the base class serialize function via base_object */
template<class Archive> 
void MovingLeastSquaresModel::serialize(Archive & archive, 
					 const unsigned int version)
{
  // serialize the base class data, then my members
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & sd;
  archive & bs;
  archive & coeffs;
  archive & continuity;
}

#endif 

#endif
