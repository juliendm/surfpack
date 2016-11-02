/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURF_POINT_H__
#define __SURF_POINT_H__

#include "surfpack_system_headers.h"
#include "SurfpackMatrix.h"

/// Holds a data point in a space of arbitrary dimension.  A SurfPoint object
/// contains an n-tuple representing the location of the point in the space, 
/// and list of zero or more response values for that point.  Includes methods
/// for querying location and responses, adding new response values, and for
/// stream I/O.
class SurfPoint {

private:

/// Nested exception class used when an attempt is made to create a SurfPoint 
/// with 0 dimensions.
class null_point : public std::runtime_error
{
public:
  null_point(const std::string& msg = 
    "Surfpack Error: attempt to make SurfPoint with 0 dimensions.") 
    : std::runtime_error(msg) {}
};

/// Nested exception for insufficient gradient data
class bad_gradient_size : public std::runtime_error
{
public:
  bad_gradient_size(const std::string& msg = 
    "Surfpack Error: SurfPoint gradient data required for all functions") 
    : std::runtime_error(msg) {}
};

/// Nested exception for insufficient Hessian data
class bad_hessian_size : public std::runtime_error
{
public:
  bad_hessian_size(const std::string& msg = 
    "Surfpack Error: SurfPoint Hessian data required for all functions") 
    : std::runtime_error(msg) {}
};
  

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:

  /// Initialize without any response values
  SurfPoint(const std::vector<double>& x);

  /// Initialize with one response value
  SurfPoint(const std::vector<double>& x, double f0);

  // The following two interfaces imply that the client has to creat
  // an extra copy of these data; could instead pass a pointer and let
  // the constructor copy it.  Or better, move Surfpack to Teuchos
  // matrices.

  /// Initialize with one response value and corresponding gradient
  SurfPoint(const std::vector<double>& x, double f0, 
	    const std::vector<double>& gradient0);

  /// Initialize with one response value and corresponding gradient and Hessian
  SurfPoint(const std::vector<double>& x, double f0, 
	    const std::vector<double>& gradient0, 
	    const SurfpackMatrix<double>& hessian0);
  
  /// Initialize with zero or more response values
  SurfPoint(const std::vector<double>& x, const std::vector<double>& f);
  
  /// Read point from istream in binary format
  SurfPoint(std::istream& is, unsigned xsize, unsigned fsize, 
	    unsigned grad_size = 0, unsigned hess_size = 0);

  /// Read point from string in text format
  SurfPoint(const std::string& single_line, unsigned xsize, unsigned fsize,
	    unsigned grad_size = 0, unsigned hess_size = 0, 
	    unsigned skip_columns = 0);

  /// Copy constructor performs a deep copy
  SurfPoint(const SurfPoint& other);

  /// Default constructor creates a one dimensional point at the origin 
  SurfPoint();

  ~SurfPoint();

private:  
  /// Initialization used by all regular constructors.  Ensures that point has
  /// at least one dimension.
  void init();


// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:

  /// Assign 'other' to 'this' unless they are already equal 
  SurfPoint& operator=(const SurfPoint& other);

  /// Tests for deep equality
  bool operator==(const SurfPoint& other) const;
  
  /// Tests for deep inequality
  bool operator!=(const SurfPoint& other) const;

  /// Return the x value along the (xindex)- th dimension (x[xindex])
  double operator[](unsigned xindex) const;
  
  /// Function object for use with pairs of SurfPoint objects (particularly in
  /// a SurfData object).  SurfPoint s1 is "less than" s2 if it has fewer
  /// dimensions.  SurfPoint s1 is also less than s2 if, for some dimension i,
  /// s1[i] < s2[i] AND s1[j] == s2[j] for all j, 0 <= j < i.  Note that the
  /// SurfPoint's response values have no bearing on the results for this 
  /// comparison.  Since the response values DO affect the results of 
  /// SurfPoint::operator==, it is NOT necessarily the case that s1 == s2 and
  /// (!SurfPointPtrLessThan(&s1,&s2) && !SurfPointPtrLessThan(&s2,&s1)) will
  /// return the same boolean value.
  class SurfPointPtrLessThan
  {
  public:
    bool operator()(const SurfPoint* sp1, const SurfPoint* sp2) const;
  };
      
// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

  /// Return dimensionality of data point
  unsigned xSize() const;

  /// Return number of response variables
  unsigned fSize() const;

  /// Return number of response variables having gradients (should be
  /// 0 or fSize())
  unsigned fGradientsSize() const;

  /// Return number of response variables having Hessians (should be 0
  /// or fSize())
  unsigned fHessiansSize() const;

  /// Return point in the domain
  const std::vector<double>& X() const;

  /// Return response value at responseIndex
  double F(unsigned responseIndex = 0) const;

  /// Return gradient vector for response value at responseIndex
  const std::vector<double>& fGradient(unsigned responseIndex = 0) const;

  /// Return response value at responseIndex
  const SurfpackMatrix<double>& fHessian(unsigned responseIndex = 0) const;


// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Append a new response variable, returning index of the inserted value
  unsigned addResponse(double val = 0); 

  /// Set an existing response variable to a new value
  void F(unsigned responseIndex, double responseValue);
 
  /// Change the value of one of the dimensions of the point
  void setX(unsigned index, double value);
 
  /// Change the dimensionality of the point
  void resize(unsigned new_size);

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write location and responses of this point to stream in binary format 
  void writeBinary(std::ostream& os) const;

  /// Write location and responses of this point to stream in text format 
  void writeText(std::ostream& os) const;

  /// Read location and responses of this point from stream in binary format 
  void readBinary(std::istream& is);
  
  /// Read location and responses of this point from stream in text format 
  void readText(const std::string& single_line, unsigned skip_columns = 0);

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected:
  /// The point in the domain.  The size of x is the dimensionality of the 
  /// space. 
  std::vector<double> x;          

  // For now, with responses, maintaining convention of a vector to
  // index responses f1, f2, ...  May later adopt continguous memory
  // data structures

  /// Zero or more response values at x (i.e., f1(x), f2(x) ... )
  std::vector<double> f;      

  /// Zero or more gradient vectors at x; indexed fGradient[k][i] =
  /// df_k/dx_i, i.e., k indexes f and i indexes x
  std::vector<std::vector<double> > fGradients;

  /// Zero or more Hessian matrices at x; 
  std::vector<SurfpackMatrix<double> > fHessians;

  
// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________
protected:

/// Provides range checking on the response values.  Throws an exception if an
/// index is requested that does not exist.
void checkRange(const std::string& header, unsigned index) const;

#ifdef __TESTING_MODE__
  friend class SurfPointTest;
  friend class SurfDataTest;
  //friend class SurfScalerTest;
#endif

private:
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class SurfPoint data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif
};

/// Write point to an output stream in text format
std::ostream& operator<<(std::ostream& os, const SurfPoint& sp); 

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SurfPoint::serialize(Archive & archive, 
			  const unsigned int version)
{  
  archive & x;          
  archive & f;      
  archive & fGradients;
  archive & fHessians;
}
#endif

#endif
