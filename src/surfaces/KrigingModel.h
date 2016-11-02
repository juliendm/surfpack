/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_H__
#define __KRIGING_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

#include "nkm/NKM_KrigingModel.hpp"

class ScaledSurfData;

/// A thin wrapper around a NewKrigingModel
class KrigingModel : public SurfpackModel
{

public:

  KrigingModel(const SurfData& sd, const ParamMap& args);
  ~KrigingModel();
  virtual double variance(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual MtxDbl hessian(const VecDbl& x) const;
  virtual std::string asString() const;

protected:

  MtxDbl getMatrix(const ScaledSurfData& ssd, const VecDbl& correlations);
  virtual double evaluate(const VecDbl& x) const;

friend class KrigingModelTest;

private:
  /// default constructor used when reading from archive file
  KrigingModel(): nkmKrigingModel(NULL) { /* empty ctor */}

  /// disallow copy construction as not implemented
  KrigingModel(const KrigingModel& other);

  /// disallow assignment as not implemented
  KrigingModel& operator=(const KrigingModel& other);

  /// helper to covert surfpack::SurfData to nkm::SurfData
  void surfdata_to_nkm_surfdata(const SurfData& sd, nkm::SurfData& nkm_sd);

  // use class data to keep in scope for wrapped model
  //nkm::KrigingModel* nkmKrigingModel;
  nkm::KrigingModel* nkmKrigingModel;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif
};


///////////////////////////////////////////////////////////
///   Kriging Model Factory	
///////////////////////////////////////////////////////////

class KrigingModelFactory : public SurfpackModelFactory 
{

public:
  KrigingModelFactory();
  KrigingModelFactory(const ParamMap& args);
  /// Override since Kriging does allow constraints
  virtual bool supports_constraints();

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  /// For Kriging, sufficient data is assessed by the NKM submodel
  virtual void sufficient_data(const SurfData& sd);
};

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
/** Serializer for the dervied Model data, e.g., basis and coefficients.
    Must call the base class serialize function via base_object */
template<class Archive> 
void KrigingModel::serialize(Archive & archive, 
			     const unsigned int version)
{
  // serialize the base class data, then my members
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & nkmKrigingModel;
}
#endif

#endif
