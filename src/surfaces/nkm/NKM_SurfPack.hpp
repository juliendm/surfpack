#ifndef __SURFPACK_HPP__
#define __SURFPACK_HPP__

#include "NKM_SurfMat.hpp"
//not sure which of these includes and usings are necessary achieved
#include <cstdlib>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <cstring>
#include <csignal>
#include <set>
#include <string>
#include <vector>

namespace nkm {

enum { SILENT_OUTPUT, QUIET_OUTPUT, NORMAL_OUTPUT, VERBOSE_OUTPUT,
       DEBUG_OUTPUT };

int if_close_enough(double a, double b);

int nchoosek(int n, int k);

int num_multi_dim_poly_coef(int nvarsr, int ndeg);

/** generate a "poly" representation of a multidimensional polynomial, each 
    column is a multidimensional monomial, for example x0^0*x1^2*x2^0*x3^1*x4^0
    is represented by [0 2 0 1 0]^T. If ndeg is negative then it produces all 
    multidimensional monomials of exactly degree |ndeg|.  If ndeg is nonnegative
    it produces all multidimensional monomials of degree up to ndeg.  0th degree
    then all 1st degree (in gradient order), then all 2nd degree (in column 
    major lower triangular part of the hessian order), then all third degree, 
    then all 4th degree, etc. Because the derivative orders in a taylor series
    approximation of a function are also its polynomial powers, this "poly" 
    representation is also well suited to representing derivative orders
*/
MtxInt& multi_dim_poly_power(MtxInt& poly, int nvarsr, int ndeg, int istart=0, int jstart=0, int iffirst=1);

/** generate a "poly" representation of main effects (no cross terms) 
    multidimensional polynomial */
MtxInt& main_effects_poly_power(MtxInt& poly, int nvarsr, int ndeg);

/** convert a "poly" representation of a polynomial to a "flypoly" (on the fly) 
    representation of a polynomial where each column is a multidimensional 
    monomial.  For example the multidimensional monomial 
    x0^0*x1^2*x2^0*x3^1*x4^0 
    has the "poly" representation [0 2 0 1 0]^T which could be evaluated as 
        coef*pow(x0,0)*pow(x1,2)*pow(x2,0)*pow(x3,1)*pow(x4,0)
    has the "flypoly" representation [3 1 1 3]^T which gets evaluated as
        coef*x1*x1*x3, the leading 3 in this flypoly example is the total 
	order or equivalently the number of multiplications needed.  
    The	"flypoly" representation has the advantage of more efficient 
    evaluation of a polynomial (when a set of points is represented as a 
    column major ordered matrix where each point is a column).  The "poly"
    representation is convenient for representing derivatives and taking 
    analytical derivatives of polynomials, the "poly" representation is used
    for "permanent" storage and "flypoly" is used for fast "on the fly" 
    purposes */
MtxInt& poly_to_flypoly(MtxInt& flypoly, const MtxInt& poly, int maxorder);

/** compute the "flypoly" representation of the polynomial that is the 
    der(:,ider) [which is a "poly" representation of a mixed multidimensional 
    derivative] of the "poly" represenation of an original polynomial */
void poly_der_to_flypoly(MtxInt& flypoly, MtxDbl& flycoef, const MtxInt& poly, 
			 const MtxInt& der, int ider, int maxorder);

/** evaluate the "flypoly" representation of polynomial at a set of points xr
    and store the result in y.  flypoly has npoly columns, coef has npoly rows
    and 1 column, xr has nvars rows and npts columns, y will have 1 row and 
    npts columns on exit */
MtxDbl& evaluate_flypoly(MtxDbl& y, const MtxInt& flypoly, const MtxDbl& coef, const MtxDbl& xr); 

/** convert the "poly" represenation of polynomial to it's "flypoly" 
    representaion by calling MtxInt& poly_to_flypoly(MtxInt& flypoly, 
    const MtxInt& poly, int maxorder); then call MtxDbl& evaluate_flypoly(
    MtxDbl& y, const MtxInt& flypoly, const MtxDbl& coef, const MtxDbl& xr); 
    to evaluate the polynomial at a set of points xr and store the result 
    in y */
MtxDbl& evaluate_poly(MtxDbl& y, MtxInt& flypoly, const MtxInt& poly, const MtxDbl& coef, const MtxDbl& xr);

/** evaluate the matrix of (multidimensional monomial, jointly a polynomial) 
    basis functions g at set of points xr from the "flypoly" representation 
    of a polynomial (where each column of xr is a point) */
MtxDbl& evaluate_flypoly_basis(MtxDbl& g, const MtxInt& flypoly, 
			       const MtxDbl& xr);

/** generate the "flypoly" represenation from the "poly" representation of a 
    polynomial and then call 
    MtxDbl& evaluate_flypoly_basis(MtxDbl& g, const MtxInt& flypoly, 
                                   const MtxDbl& xr) 
    to evaluate the matrix of basis functions g at set of points xr (where 
    each column of xr is a point)*/
MtxDbl& evaluate_poly_basis(MtxDbl& g, MtxInt& flypoly, const MtxInt& poly, 
			    const MtxDbl& xr);


/** analytically take the "der" set of mixed multidimensional derivatives 
    (stored in a "poly" style representation) of the "poly" reprensentation
    of an original polynomial and then evaluate this set of derivatives at a
    set of points xr and store the result in dy.  poly has nvars rows and
    npoly columns, coef has npoly rows and 1 column, der has nvars rows and 
    nder columns, xr has nvars rows and npts columns.  dy will have nder rows
    and npts columns on exit. flypoly and flycoef are workspace variables */
MtxDbl& evaluate_poly_der(MtxDbl& dy, MtxInt& flypoly, MtxDbl& flycoef, 
			  const MtxInt& poly, const MtxInt& der, 
			  const MtxDbl& coef, const MtxDbl& xr);

/** evaulate a matrix dg containing the "der" set mixed multidimensional 
    derivatives of the "poly" representation of an original polynomial at 
    at set of points xr (where each column of xr is a point), flypoly and 
    flycoef are work space variables */
MtxDbl& evaluate_poly_der_basis(MtxDbl& dg, MtxInt& flypoly, MtxDbl& flycoef, 
				const MtxInt& poly, const MtxInt& der, 
				const MtxDbl& xr);

/// generate a rotation matrix from a set of Euler angles 
MtxDbl& gen_rot_mat(MtxDbl& Rot, const MtxDbl& EulAng, int nvarsr);

/// generate a random set of Euler angles and convert them to a rotation matrix
MtxDbl& gen_rand_rot_mat(MtxDbl& rot,int nvarsr);

///generates 2*nvarsr random samples between 0 and 1, the sample design is stored in a matrix with nvarsr rows and 2*nvars columns (one point per column), the sample design is binning optimal with the bins being chosen as the end points of a randomly rotated set of axes (so the BINS, but not the points, are maximin spaced) the design is NOT a latin hypercube and it is not symmetric, opposite octants are stored in sequential columns of the xr matrix.
MtxDbl& gen_rand_axis_bin_opt_samples_0to1(MtxDbl& xr, int nvarsr); 

//xr is for XReal, nvarsr is for Number of VARiables Real
inline MtxDbl& rotate_xr(MtxDbl& xr_rot, const MtxDbl& rot_or_eul_ang, const MtxDbl& xr) {
    int nvarsr=xr.getNCols();
    if((rot_or_eul_ang.getNRows()==nvarsr)&&
       (rot_or_eul_ang.getNCols()==nvarsr)) {
      //rot_or_eul_ang is a rotation matrix
      matrix_mult(xr_rot,rot_or_eul_ang,xr,0.0,1.0,'T','N');
    }
    else if((rot_or_eul_ang.getNRows()==((nvarsr*(nvarsr-1))/2))&&
	    (rot_or_eul_ang.getNCols()==1)) {
      //rot_or_eul_ang is a vector containing the Euler Angles for the rotation matrix
      MtxDbl rot; 
      gen_rot_mat(rot,rot_or_eul_ang,nvarsr);
      matrix_mult(xr_rot,rot,xr,0.0,1.0,'T','N');
    }
    else{
      std::printf("Error in rotate_xr(MtxDbl& xr_rot,MtxDbl& rot_or_eul_ang,MtxDbl& xr): rot_or_eul_ang has the wrong size!!!\n");
      assert(0);
    }
    return xr_rot;
}


inline double dsign(double x){return static_cast<double>((x>0.0)-(x<0.0));};
inline int isign(int x){return static_cast<int>((x>0)-(x<0));};


template<typename T>
std::vector<T>& toVec(std::vector<T>& result, const std::string& s)
{
  std::istringstream is(s);
  result.clear();
  if (s == "") 
    return result;
  T temp;
  do {
    is >> temp;
    result.push_back(temp);
  } while (!is.eof());
  
  return result;
}


template<typename T>
std::string toString(const T arg)
{
  std::ostringstream os;
  os << arg;
  return os.str();
}

template<typename T>
std::string fromVec(const std::vector<T>& vec)
{
  std::ostringstream os;
  for (typename std::vector<T>::const_iterator itr = vec.begin();
	itr != vec.end(); ++itr) {
    if (itr != vec.begin()) os << " ";
    os << *itr;
  }
  return os.str();
}

namespace surfpack {
  ///should have variable # of sig figs (precision) control for output

  /// Precision of output for double precision numbers
  const unsigned output_precision = 16;

  /// Length of the field for double-precision number stream output
  const unsigned field_width = output_precision + 6;


  /// Thrown when an attempt to open a file for reading or writing fails
  class file_open_failure: public std::runtime_error
  {
  public:
    file_open_failure(const std::string& filename = "") 
      : std::runtime_error("File " + filename + " could not be opened.") {}
  };
    
  /// Thrown when end-of-file is reached unexpectedly, when an unrecognized or
  /// unacceptable file extension is encountered, or when a file contains
  /// unexpected or illegally formatted contents.
  class io_exception: public std::runtime_error
  {
  public:
    io_exception(const std::string& msg = "") : std::runtime_error(msg) {}
  };

  /// Throw an exception if end-of-file has been reached 
  void checkForEOF(std::istream& is);

  /// specified by parameter extension
  bool hasExtension(const std::string& filename, const std::string extension);
}

} // end namespace nkm

#endif
