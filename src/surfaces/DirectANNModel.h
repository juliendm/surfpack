/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __DIRECT_ANN_MODEL_H__
#define __DIRECT_ANN_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"


/// The basis set is the first layer in the neural network.  It
/// applies random weights and bias to the (previously scaled) inputs
/// and a tanh activation function.  gamma(x) = tanh( A0 x + theta0 )
class DirectANNBasisSet
{
public:
  /// random weights and bias [ A0 | theta0 ] for input layer
  MtxDbl weights;

  /// default constructor used when reading from archive file
  DirectANNBasisSet() { /* empty ctor */ }

  /// standard constructor accepting weights [ A0 | theta0 ] to apply
  /// to the input layer
  DirectANNBasisSet(const MtxDbl& weights_in);
  
  /// evaluate the index-th basis function at the point x
  double eval(unsigned index, const VecDbl& x) const;

  /// evaluate the index-th basis function derivative at the point x
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;

  /// compute the contribution due to the index-th basis function at the point x
  double nodeSum(unsigned index, const VecDbl& x) const;

  /// write the basis set as a string
  std::string asString() const;

private:

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private serialize function
  friend class boost::serialization::access;
  /// serializer for derived class data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};

/** The model consists of a basis set together with coefficients
    (typically estimated from data via least squares) and a tanh
    activation function.  The final model is 
    tanh( A1 tanh( A0 x + theta0 ) + theta1 ) */
class DirectANNModel : public SurfpackModel
{

public:

  DirectANNModel(const DirectANNBasisSet& bs_in, const VecDbl& coeffs_in);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;

protected:

  /// default constructor used when reading from archive file
  DirectANNModel() { /* empty ctor */ }

  /// evaluate the model at the point x
  virtual double evaluate(const VecDbl& x) const;

  /// basis set mapping the input layer to the hidden layer
  DirectANNBasisSet bs;

  /// coefficients and bias [ A1 | theta1 ] for hidden layer
  VecDbl coeffs;


private:

  /// disallow copy construction as not implemented
  DirectANNModel(const DirectANNModel& other);

  /// disallow assignment as not implemented
  DirectANNModel& operator=(const DirectANNModel& other);

friend class DirectANNModelTest;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};

///////////////////////////////////////////////////////////
///   Direct ANN Model Factory	
///////////////////////////////////////////////////////////


/** Factory to generate neural network models.  If user requests a
    specific number of nodes, use that number, up to a maximum of
    data.size()-1.  If no specific request, use data.size()-1, with a
    cap at 100 and use OMP to downselect. */
class DirectANNModelFactory : public SurfpackModelFactory 
{

public:
  DirectANNModelFactory();
  DirectANNModelFactory(const ParamMap& args);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  /// generate a random matrix over (-range/2, range/2)
  MtxDbl randomMatrix(unsigned nrows, unsigned ncols);

  /// number of user-requested hidden layer nodes
  unsigned maxNodes;
  /// range of the random coefficients for inner layer (default 2.0 --> (-1, 1))
  double range;
  /// unused
  unsigned samples;
  /// random seed for cross validation-based basis selection
  unsigned randomSeed;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  /// allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

template< class Archive >
void DirectANNBasisSet::serialize(Archive & archive, 
				  const unsigned int version)
{
  archive & weights;
}

template< class Archive >
void DirectANNModel::serialize(Archive & archive,
			       const unsigned int version)
{
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & bs;
  archive & coeffs;
}

#endif

#endif
