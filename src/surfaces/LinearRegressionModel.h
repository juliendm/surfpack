/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __LINEAR_REGRESSION_MODEL_H__
#define __LINEAR_REGRESSION_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

class SurfPoint;
class ScaledSurfData;

class LRMBasisSet
{
public:

  VecVecUns bases;
  double eval(unsigned index, const VecDbl& x) const;
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;
  std::string asString() const;
  void add(const std::string& s_basis);
  unsigned size() const { return bases.size();}

private:

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for linear regression bases
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


class LinearRegressionModel : public SurfpackModel
{
public:

  /// standard constructor from a basis set
  LinearRegressionModel(const unsigned dims, const LRMBasisSet& bs_in, 
			const VecDbl& coeffs_in);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;

protected:

  /// default constructor used when reading from archive file 
  LinearRegressionModel() { /* empty ctor */ }

  virtual double evaluate(const VecDbl& x) const;
  LRMBasisSet bs;
  VecDbl coeffs;

friend class LinearRegressionModelTest;

private:

  /// disallow copy construction as not implemented
  LinearRegressionModel(const LinearRegressionModel& other);

  /// disallow assignment as not implemented
  LinearRegressionModel& operator=(const LinearRegressionModel& other);

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


struct Term {
  bool color;
  VecUns vars;
  Term(const VecUns& vars_in) : color(false), vars(vars_in) {}
};
  

///////////////////////////////////////////////////////////
///	Linear Regression Model Factory	
///////////////////////////////////////////////////////////

class LinearRegressionModelFactory : public SurfpackModelFactory 
{

public:

  LinearRegressionModelFactory();
  LinearRegressionModelFactory(const ParamMap& args);
  virtual unsigned minPointsRequired();
  virtual unsigned recommendedNumPoints();
  /// LRM does allow constraints
  virtual bool supports_constraints();
  VecDbl lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd);
  static LRMBasisSet CreateLRM(unsigned order, unsigned dims);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  /// Sufficient data is based on points plus constraint data
  virtual void sufficient_data(const SurfData& sd);
  unsigned order;
  MtxDbl eqConLHS;
  VecDbl eqConRHS;

private:

 /// convenience function to create constraint linear system in the factory
  void setEqualityConstraints(const SurfPoint& sp);

};


// -----
// Definitions and export of serialization functions
// -----

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

template<class Archive> 
void LRMBasisSet::serialize(Archive & archive, const unsigned int version)
{
  archive & bases;
}

/** Serializer for the derived Model data, e.g., basis and coefficients.
    Must call the base class serialize function via base_object */
template<class Archive> 
void LinearRegressionModel::serialize(Archive & archive, 
				      const unsigned int version)
{
  // serialize the base class data, then my members
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & bs;
  archive & coeffs;
}


#endif


#endif  // __LINEAR_REGRESSION_MODEL_H__
