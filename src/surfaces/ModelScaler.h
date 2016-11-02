/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MODEL_SCALER_H__
#define __MODEL_SCALER_H__

#include "surfpack_system_headers.h"

class SurfData;

class ModelScaler {

public:

  virtual const VecDbl& scale(const VecDbl& unscaled_x) const = 0;
  virtual double descale(double scaled_response) const = 0;
  virtual double scaleResponse(double unscaled_response) const = 0;
  virtual std::string asString() { return ""; }
  ModelScaler() {}
  virtual ~ModelScaler() {}
  virtual ModelScaler* clone() const = 0;

private:

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private serializer
  friend class boost::serialization::access;
  /// serializer for base class Model data
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version);
#endif

};

class NonScaler : public ModelScaler {

public:

  virtual const VecDbl& scale(const VecDbl& unscaled_x) const;
  virtual double descale(double scaled_response) const ;
  virtual double scaleResponse(double unscaled_response) const ;
  virtual std::string asString();
  static ModelScaler* Create(const SurfData& data);
  virtual ModelScaler* clone() const;

private:

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private serializer
  friend class boost::serialization::access;
  /// serializer for base class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


class NormalizingScaler : public ModelScaler {

public:
  struct Scaler {
    double offset;
    double scaleFactor;
    Scaler(double o, double s) : offset(o), scaleFactor(s) {}
    Scaler() {}
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    // allow serializers access to private data
    friend class boost::serialization::access;
    /// serializer for base class Model data
    template<class Archive> 
    void serialize(Archive & archive, const unsigned int version);
#endif
  };

  virtual const VecDbl& scale(const VecDbl& unscaled_x) const;
  virtual double descale(double scaled_response) const;
  virtual double scaleResponse(double unscaled_response) const ;
  virtual std::string asString();
  virtual VecDbl getScalerOffsets() const;
  virtual VecDbl getScalerScaleFactors() const;
  virtual double getDescalerOffset() const;
  virtual double getDescalerScaleFactor() const;
  NormalizingScaler(const std::vector<Scaler>& s, const Scaler& d) 
    : scalers(s), descaler(d), result(s.size()) {}
  ~NormalizingScaler() {}
  // constructor to normalize each var/resp to [ 0, 1 ]
  static ModelScaler* Create(const SurfData& data);
  // constructor to normalize each var/resp to [ -norm_factor, norm_factor ]
  static ModelScaler* Create(const SurfData& data, double norm_factor);
  virtual ModelScaler* clone() const;
  
protected:

  std::vector<Scaler> scalers;
  Scaler descaler;
  mutable VecDbl result;
  friend class ModelScalerTest;

private:

  NormalizingScaler() { }

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for base class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};

class ScaledSurfData {

public:

  ScaledSurfData(const ModelScaler& ms_in, const SurfData& sd_in);
  std::vector< double > getResponses() const;
  double getResponse(unsigned index) const;
  unsigned size() const;
  unsigned xSize() const;
  double operator()(unsigned pt, unsigned dim) const;
  const std::vector<double>& operator()(unsigned pt) const;
  static VecVecDbl asVecVecDbl(const ScaledSurfData& data);

protected:

  const ModelScaler& ms;
  const SurfData& sd;

};


// -----
// Definitions and export of serialization functions
// -----

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

template<class Archive>
void NormalizingScaler::Scaler::serialize(Archive & archive, 
					  const unsigned int version)
{
  archive & offset;
  archive & scaleFactor;
}


/** Serializer for the base class ModelScaler data.  This is
    called by the derived class serialize functions via base_object */
template<class Archive>
void ModelScaler::serialize(Archive & archive, const unsigned int version)
{
  // Nothing to serialize at base
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(ModelScaler)


template<class Archive> 
void NonScaler::serialize(Archive & archive, const unsigned int version)
{
  // serialize the base class data, then my members
  archive & boost::serialization::base_object<ModelScaler>(*this);
  // Nothing to serialize
}


template<class Archive> 
void NormalizingScaler::serialize(Archive & archive, const unsigned int version)
{
  // serialize the base class data, then my members 
  archive & boost::serialization::base_object<ModelScaler>(*this);
  archive & scalers;
  archive & descaler;
  archive & result;
}

#endif

#endif  // __MODEL_SCALER_H__

