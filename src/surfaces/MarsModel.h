/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MARS_MODEL_H__
#define __MARS_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

typedef float real;


class MarsModel : public SurfpackModel
{

public:

  /// standard constructor from a factory
  MarsModel(const unsigned dims, real* fm_in, int fmsize, int* im_in, 
    int imsize, int interp);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;

protected:

  /// default constructor used when reading from archive file 
  MarsModel() { /* empty ctor */ }

  virtual double evaluate(const VecDbl& x) const;

  std::vector<real> fm;
  std::vector<int> im;
  int interpolation;

friend class MarsModelTest;

private:

  /// disallow copy construction as not implemented
  MarsModel(const MarsModel& other);

  /// disallow assignment as not implemented
  MarsModel& operator=(const MarsModel& other);

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


///////////////////////////////////////////////////////////
///	Linear Regression Model Factory	
///////////////////////////////////////////////////////////

class MarsModelFactory : public SurfpackModelFactory 
{

public:

  MarsModelFactory();
  MarsModelFactory(const ParamMap& args);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  real* xMatrix;
  real* fm;
  int* im;
  int n;
  int np;
  int max_bases; 
  int max_interactions;
  int interpolation;
};


// -----
// Definitions and export of serialization functions
// -----

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

/** Serializer for the derived Model data, e.g., basis and coefficients.
    Must call the base class serialize function via base_object */
template<class Archive>
void MarsModel::serialize(Archive & archive, const unsigned int version)
{
  // serialize the base class data, then my members
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & fm;
  archive & im;
  archive & interpolation;
}

#endif


#endif  // __MARS_MODEL_H__


