/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_H__
#define __SURFPACK_H__

#include "surfpack_system_headers.h"

#include "MersenneTwister.h"
#include "SurfpackMatrix.h"

//class AbstractSurfDataIterator;
class SurfData;

enum DifferenceType {
  DT_ABSOLUTE,
  DT_SQUARED,
  DT_SCALED
};

enum MetricType {
  MT_RELATIVE_MAXIMUM,
  MT_RELATIVE_AVERAGE,
  MT_MINIMUM,
  MT_MAXIMUM,
  MT_SUM,
  MT_MEAN,
  MT_ROOT_MEAN
};

namespace surfpack {

enum { SILENT_OUTPUT, QUIET_OUTPUT, NORMAL_OUTPUT, VERBOSE_OUTPUT,
       DEBUG_OUTPUT };

// _____________________________________________________________________________
// Debugging Output Strategy 
// _____________________________________________________________________________

class DbgStream {
public:
  mutable int level;
  DbgStream() : level(0) {}
  ~DbgStream() {}
  const DbgStream& operator()(int level_in) const {
    level = level_in;
    return *this;
  }
  template<typename T> const DbgStream& operator<<(const T& item) const {
    if (level) {
      std::cout << item;
    }
    return *this;
  }
};
const DbgStream& dbg(int level_in);
// _____________________________________________________________________________
// Mersenne Twister Random Number Generator 
// _____________________________________________________________________________
                                                                                
class MyRandomNumberGenerator : std::unary_function<int,int>
{
public:
  MyRandomNumberGenerator() {}
  MTRand mtrand;
  // operator for std::random_shuffle; int in [0, n-1]
  int operator()(int n)
  {
    return mtrand.randInt(n-1);
  }
  void seed(int seeder)
  {
    mtrand.seed(seeder);
  }
  /// double in [0,1]
  double rand()
  {
    return mtrand.rand();
  }
  /// double in [0,1)
  double randExc()
  {
    return mtrand.randExc();
  }
  /// int in [0,n]
  int randInt(int n)
  {
    return mtrand.randInt(n);
  }

};
                                                                                
MyRandomNumberGenerator& shared_rng();

// _____________________________________________________________________________
// Block partitioning helper methods 
// _____________________________________________________________________________

unsigned block_low(unsigned id, unsigned p, unsigned n);
unsigned block_high(unsigned id, unsigned p, unsigned n);
unsigned block_size(unsigned id, unsigned p, unsigned n);
unsigned block_owner(unsigned j, unsigned p, unsigned n);

// _____________________________________________________________________________
// Constants 
// _____________________________________________________________________________
  // Precision of output for double precision numbers
  const unsigned output_precision = 6;

  // Length of the field for double-precision number stream output
  const unsigned field_width = output_precision + 9;

// _____________________________________________________________________________
// Nested Types 
// _____________________________________________________________________________

  /// For use in comparing the actual value of some function with the estimate
  /// given by a Surface approximation
  struct ErrorStruct {
    double observed;
    double estimated;
  };

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
  
// ____________________________________________________________________________
// MISC functions
// ____________________________________________________________________________

// windows doesn't have a native atanh function, long term we may want to switch to boost but for now we implement it ourselves.
  inline double atanh(double x)
  {
#if defined(_WIN32) || defined(_WIN64)
    return 0.5*(std::log(1.0+x) - std::log(1.0-x));
#else
    return ::atanh(x);
#endif
  };


// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________
  
  /// Write the value of contents to the file specified by filename.  Throw an
  /// exception if the file cannot be opened.
  void writeFile(std::string filename, std::string contents);

  /// Write the parameter header, followed by the matrix mat (the dimensions of
  /// which are specified by parameters rows and columns) to the parameter os.
  /// If c_style is true, the memory layout is assumed to follow the C
  /// convention (if mat points to an m by n matrix, the first m values are
  /// interpreted as the first row).  Otherwise, the layout is assumed to 
  /// follow the Fortran convention (the first n values are interpreted as the 
  /// first column).
  void writeMatrix(const std::string header, double* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);

  /// Write the parameter header, followed by the matrix mat (the dimensions of
  /// which are specified by parameters rows and columns) to the parameter os.
  /// If c_style is true, the memory layout is assumed to follow the C
  /// convention (if mat points to an m by n matrix, the first m values are
  /// interpreted as the first row).  Otherwise, the layout is assumed to 
  /// follow the Fortran convention (the first n values are interpreted as the 
  /// first column).
  void writeMatrix(const std::string header, unsigned* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);

  /// Write the contents of a matrix to a file specified by parameter filename.
  /// Open the file and call another version of writeMatrix.
  void writeMatrix(const std::string filename, double* mat, unsigned rows, 
    unsigned columns, bool c_style = false);

  /// Write the contents of a matrix to a file specified by parameter filename.
  /// Open the file and call another version of writeMatrix.
  void writeMatrix(const std::string filename, unsigned* mat, unsigned rows, 
    unsigned columns, bool c_style = false);

  /// Write the parameter header followed by the values in the vector
  void printVector(const std::string header, VecDbl& vec, 
    std::ostream& os = std::cout);

  /// Return true if the file specified by parameter file name has the extension
  /// specified by parameter extension
  bool hasExtension(const std::string& filename, const std::string extension);

  /// Return true for binary model filename, false for text
  bool isBinaryModelFilename(const std::string& filename);

  /// Throw an exception if end-of-file has been reached 
  void checkForEOF(std::istream& is);

  /// Open the file specified by filename and return the type of Surface that 
  /// is written there.  Throw an exception if the file cannot be opened, or if
  /// the file extension is anything other than .txt or .srf. 
  const std::string surfaceName(const std::string filename);

  /// Return the next item in the file as a string.  If the file is opened in
  /// binary mode, first read an integer that specifies the number of 
  /// characters in the string, then read the string. 
  const std::string readName(std::istream& is, bool binary);

  /// Round values that are close to integers to integers
  void approximateByIntegers(VecDbl& vals, double epsilon = 1.e-6);

// ____________________________________________________________________________
// Vector helper methods 
// ____________________________________________________________________________

  /// Return the sum of the vector of values
  double sum_vector(VecDbl& vals);

  /// Return the arithmetic mean (average) of the values in vector vals
  double mean(const VecDbl& vals);

  /// Return the sample variance of the values in vals
  double sample_var(VecDbl& vals);
  
  /// Return the sample standard deviation of the values in vals
  double sample_sd(VecDbl& vals);

  /// Return the sum of squared deviations from the mean
  double sum_squared_deviations(VecDbl& vals);

  /// Return the sum of absolute deviations from the mean
  double sum_absolute_deviations(VecDbl& vals);

  /// Return absolute, squared, or relative differences of second and third
  /// parameters through the first parameter
  void differences(VecDbl& results, VecDbl& observed,
    VecDbl& predicted, enum DifferenceType dp = DT_ABSOLUTE);
  
  /// Return the euclidean distance between pt1 and pt2.  Throw an exception if
  /// the dimensionality of the two vectors does not match.
  double euclideanDistance(const VecDbl& pt1, 
    const VecDbl& pt2);

  /// Store the vector difference between pt1 and pt2 in the paramter diff.
  /// Throw an exception if the dimensionality of the points does not match.
  void vectorDifference(VecDbl& diff, const VecDbl& pt1, const VecDbl& pt2);
// ____________________________________________________________________________
// Functions for common linear algebra tasks 
// ____________________________________________________________________________
  /// Least squares solve of system Ax = b
  void linearSystemLeastSquares(MtxDbl& A, VecDbl& x, VecDbl b);

  /// Least squares solve os system Ax = c, subject to Bx = d
  void leastSquaresWithEqualityConstraints(MtxDbl& A, 
    VecDbl& x, VecDbl& c,
    MtxDbl& B, VecDbl& d);

  /// Calls dgetrf followed by dgetri
  MtxDbl& inverse(MtxDbl& matrix);

  /// Calls dgetrf to compute LU Decomposition
  MtxDbl& LUFact(MtxDbl& matrix, 
    std::vector<int>& ipvt);

  /// Calls dgetri to compute matrix inverse, after prior call to dgetrf
  MtxDbl& inverseAfterLUFact(MtxDbl& matrix, std::vector<int>& ipvt);

  /// Note: These matrix functions would not fit easily in SurfpackMatrix.h
  /// because the fortran math functions are not templated
  /// matrix-vector mutltiplication
  VecDbl& matrixVectorMult(VecDbl& result,
    MtxDbl& matrix, VecDbl& the_vector,
    char trans = 'N');

  /// matrix-matrix multiplication
  MtxDbl& matrixMatrixMult(MtxDbl& result, MtxDbl& matrixA, MtxDbl& matrixB,
    char transA = 'N', char transB = 'N');

  /// matrix-matrix addition
  MtxDbl& matrixSum(MtxDbl& result, MtxDbl& matrixA, MtxDbl& matrixB);
  MtxDbl& matrixSubtraction(MtxDbl& result, MtxDbl& matrixA, MtxDbl& matrixB);

  /// vector-vector inner product
  double dot_product(const VecDbl& vector_a, const VecDbl& vector_b);

  /// Adds or subtracts same value to all vector elements
  VecDbl& vectorShift(VecDbl& the_vector, double shift_value);
  
  /// Returns the weighted average of two vectors: alpha*first+(1-alpha)*second
VecDbl weightedAvg(const VecDbl& first, const VecDbl& second, double alpha = 0.5);
// ____________________________________________________________________________
// Converting a string to a vector 
// ____________________________________________________________________________
template<typename T>
std::vector<T> toVec(const std::string& s)
{
  std::istringstream is(s);
  std::vector<T> result;
  if (s == "") return result;
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

void stripQuotes(std::string& str);

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________
  /// Return the value of the test function specified by parameter name at the 
  /// point specified by parameter pt
  /// \todo Change if-else construct to lookup in STL map of <name, function>
  /// pairs.
  double testFunction(const std::string name, const VecDbl& pt);

  /// Non-trivial polynomial function
  double moderatepoly(const VecDbl& pt);

  /// Tony Giunta's test function: a little wave mixed with a big wave, plus 
  /// noise
  double quasisine(const VecDbl& pt);

  /// f(x) = sigma{i=1 to n}(x_i^2 - 10*cos(2*pi*x) + 10).  With side 
  /// constraints near zero (e.g. +/-10) along each dimension, the function 
  /// appears highly-multimodal.  With larger bounds, the bowl shape becomes
  /// more dominant and the waviness is reduced to noise.
  double rastrigin(const VecDbl& pt);

  /// A multi-dimensional extension of the classic Rosenbrock test function
  double rosenbrock(const VecDbl& pt);

  /// f(x) = 3 + sigma{i=1 to n}(2*x_i)
  double simplepoly(const VecDbl& pt);

  /// Sum of the sine function along each dimension
  double sinewave(const VecDbl& pt);

  /// Sum of squares along each dimension
  double sphere(const VecDbl& pt);

  /// f(x) = sigma{i=1 to i}(x_i)
  double sumofall(const VecDbl& pt);

  /// f(x) = sigma{i=1 to n}(x_i + sin x_i)
  double xplussinex(const VecDbl& pt);

  /// Random (different queries for the same point will give different results)
  double noise(const VecDbl& pt);
} // namespace surfpack
 
#endif // __SURFPACK_H__
