/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "surfpack_LAPACK_wrappers.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::ios;
using std::numeric_limits;
using std::ofstream;
using std::ostream;
using std::setw;
using std::string;
using std::vector;


// BMA: Debugging KrigingModel serialization...
// #include "KrigingModel.h"
// #include <boost/archive/text_iarchive.hpp>
// #include <boost/archive/text_oarchive.hpp>

// template void KrigingModel::serialize<boost::archive::text_iarchive>(
//     boost::archive::text_iarchive & ar, 
//     const unsigned int file_version
// );
// template void KrigingModel::serialize<boost::archive::text_oarchive>(
//     boost::archive::text_oarchive & ar, 
//     const unsigned int file_version
// );

// _____________________________________________________________________________
// Debugging Output Strategy 
// _____________________________________________________________________________
const surfpack::DbgStream& surfpack::dbg(int level_in)
{
  static DbgStream dbg_stream;
  return dbg_stream(level_in);
}

// _____________________________________________________________________________
// Mersenne Twister Random Number Generator 
// _____________________________________________________________________________
surfpack::MyRandomNumberGenerator& surfpack::shared_rng()
{
  static MyRandomNumberGenerator mrng;
  return mrng;
}
                                                                           
// _____________________________________________________________________________
// Block partitioning helper methods 
// _____________________________________________________________________________

unsigned surfpack::block_low(unsigned id, unsigned p, unsigned n)
{
  return id*n/p;
}

unsigned surfpack::block_high(unsigned id, unsigned p, unsigned n)
{
  return block_low(id+1,p,n)-1;
}

unsigned surfpack::block_size(unsigned id, unsigned p, unsigned n)
{
  return block_high(id,p,n)-block_low(id,p,n)+1;
}

unsigned surfpack::block_owner(unsigned j, unsigned p, unsigned n)
{
  return (p*(j+1)-1)/n;
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

/// Write the value of contents to the file specified by filename.  Throw an
/// exception if the file cannot be opened.
void surfpack::writeFile(string filename, string contents)
{
  ofstream outfile(filename.c_str(), ios::out);
  if (!outfile) {
    throw surfpack::file_open_failure(filename);
  }
  outfile << contents << endl;
  outfile.close();
}


/// Write the parameter header, followed by the matrix mat (the dimensions of
/// which are specified by parameters rows and columns) to the parameter os.
/// If c_style is true, the memory layout is assumed to follow the C
/// convention (if mat points to an m by n matrix, the first m values are
/// interpreted as the first row).  Otherwise, the layout is assumed to 
/// follow the Fortran convention (the first n values are interpreted as the 
/// first column).
void surfpack::writeMatrix(const string header, double* mat, unsigned rows, 
  unsigned columns, ostream& os, bool c_style)
{
  if (header != "none" && header != "") {
    os << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(surfpack::field_width) << mat[c + r * columns];
      } else {
        os << setw(surfpack::field_width) << mat[r + c * rows];
      }
    }
    os << endl;
  }
}

/// Write the parameter header, followed by the matrix mat (the dimensions of
/// which are specified by parameters rows and columns) to the parameter os.
/// If c_style is true, the memory layout is assumed to follow the C
/// convention (if mat points to an m by n matrix, the first m values are
/// interpreted as the first row).  Otherwise, the layout is assumed to 
/// follow the Fortran convention (the first n values are interpreted as the 
/// first column).
/// \todo Templatize the method.  Priority: low.
void surfpack::writeMatrix(const string header, unsigned* mat, unsigned rows, 
  unsigned columns, ostream& os, bool c_style)
{
  if (header != "none" && header != "") {
    os << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(surfpack::field_width) << mat[c + r * columns];
      } else {
        os << setw(surfpack::field_width) << mat[r + c * rows];
      }
    }
    os << endl;
  }
}

/// Write the contents of a matrix to a file specified by parameter filename.
/// Open the file and call another version of writeMatrix.
void surfpack::writeMatrix(const string filename, double* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    throw file_open_failure(filename);
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

/// Write the contents of a matrix to a file specified by parameter filename.
/// Open the file and call another version of writeMatrix.
void surfpack::writeMatrix(const string filename, unsigned* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    throw file_open_failure(filename);
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

/// Write the parameter header followed by the values in the vector
/// \todo Use an output iterator instead.  Priority: very low.
void surfpack::printVector(const string header, vector<double>& vec,
  ostream& os)
{
  os << header << " size: " << vec.size() << endl;
  for (unsigned i = 0; i < vec.size(); i++) {
    os << i << " " << vec[i] << endl;
  }
}

/// Return true if the file specified by parameter file name has the extension
/// specified by parameter extension
bool surfpack::hasExtension(const string& filename, const string extension)
{
  return (filename.find(extension) == filename.size() - extension.size());
}

/// Validate model (surface) filename.  Return true if filename has
/// .bsps extension, false if .sps extension, otherwise throws
/// surfpack::io_exception.
bool surfpack::isBinaryModelFilename(const string& filename)
{
  bool binary;
  if (surfpack::hasExtension(filename,".bsps")) {
    binary = true;;
  } else if (surfpack::hasExtension(filename,".sps")) {
    binary = false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized model (surface) filename extension.  Use .sps or .bsps"
    );
  }
  return binary;
}


/// Throw an exception if end-of-file has been reached 
void surfpack::checkForEOF(istream& is)
{
  if (is.eof()) {
    throw surfpack::io_exception("End of file reached unexpectedly.");
  }
}

/// Open the file specified by filename and return the type of Surface that 
/// is written there.  Throw an exception if the file cannot be opened, or if
/// the file extension is anything other than .sps or .bsps. 
const string surfpack::surfaceName(const string filename)
{
  bool binary = isBinaryModelFilename(filename);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  } 
  string nameInFile = readName(infile, binary);  
  infile.close();
  return nameInFile;
} 

/// Return the next item in the file as a string.  If the file is opened in
/// binary mode, first read an integer that specifies the number of 
/// characters in the string, then read the string. 
const string surfpack::readName(istream& is, bool binary)
{
  string nameInFile;
  if (!binary) {
    getline(is,nameInFile);
    return nameInFile;
  } else {
    unsigned nameSize;
    is.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    char* surfaceType = new char[nameSize+1];
    is.read(surfaceType,nameSize);
    surfaceType[nameSize] = '\0';
    return string(surfaceType);
  }
}

/// Round values that are close to integers to integers
void surfpack::approximateByIntegers(vector<double>& vals, double epsilon)
{
  for(vector<double>::iterator iter = vals.begin(); iter != vals.end();
    ++iter) {
    double approx = static_cast<double>(static_cast<int>(*iter));
    if (fabs(*iter-approx) < epsilon) {
      *iter = approx;
    }
  }
}
// ____________________________________________________________________________
// Vector helper methods 
// ____________________________________________________________________________

/// Return the sum of the vector of values
double surfpack::sum_vector(vector<double>& vals)
{
  double sum = 0;
  for (unsigned i = 0; i < vals.size(); i++) {
    sum += vals[i];
  }
  return sum;
}

/// Return the arithmetic mean (average) of the values in vector vals
double surfpack::mean(const VecDbl& vals)
{
  //return sum_vector(vals) / vals.size();
  // Check for 'inf' and ignore those values
  unsigned excludedInfs = 0;
  double sum = 0;
  for (unsigned i = 0; i < vals.size(); i++) {
    // false when vals[i] == inf or nan
    if ( vals[i] != numeric_limits<double>::infinity()) { 
      sum += vals[i];
    } else {
      excludedInfs++;
    }
  }
  return sum / (vals.size() - excludedInfs);
}

/// Return the sample variance of the values in vals
double surfpack::sample_var(vector<double>& vals)
{
  double sse = sum_squared_deviations(vals);
  return sse / (vals.size() - 1);
}

/// Return the sample standard deviation of the values in vals
double surfpack::sample_sd(vector<double>& vals)
{
  return sqrt(surfpack::sample_var(vals));
}

/// Return the sum of squared deviations from the mean
double surfpack::sum_squared_deviations(vector<double>& vals)
{
  double sse = 0;
  double avg = surfpack::mean(vals);
  for (unsigned i = 0; i < vals.size(); i++) {
    sse += (vals[i]-avg)*(vals[i]-avg);
  }
  return sse;
}
  
/// Return the sum of absolute deviations from the mean
double surfpack::sum_absolute_deviations(vector<double>& vals)
{
  double sae = 0;
  double avg = surfpack::mean(vals);
  for (unsigned i = 0; i < vals.size(); i++) {
    sae += (vals[i]-avg);
  }
  return sae;
}
/// Return absolute, squared, or relative differences of second and third
/// parameters through the first parameter
void surfpack::differences(vector<double>& results, 
  vector<double>& observed, vector<double>& predicted, 
  enum DifferenceType dp)
{
  results.resize(observed.size());
  for (unsigned i = 0; i < observed.size(); i++) {
    results[i] = fabs(observed[i] - predicted[i]);
    switch (dp) {
      case DT_SQUARED: results[i] *= results[i]; break;
      case DT_SCALED: results[i] /= fabs(observed[i]); break;
    }
  }
}

/// Return the euclidean distance between pt1 and pt2.  Throw an exception if
/// the dimensionality of the two vectors does not match.
double surfpack::euclideanDistance(const vector<double>& pt1, 
  const VecDbl& pt2)
{
  double distance = 0.0;
  if (pt1.size() != pt2.size()) {
    throw string(
      "Cannot compute euclidean distance.  Vectors have different sizes.");
  } else {
    for (unsigned i = 0; i < pt1.size(); i++) {
      distance += (pt1[i] - pt2[i]) * (pt1[i] - pt2[i]);
    }
    distance = sqrt(distance);
  }
  return distance;
}

/// Store the vector difference between pt1 and pt2 in the paramter diff.
/// Throw an exception if the dimensionality of the points does not match.
void surfpack::vectorDifference(vector<double>& diff, const vector<double>& pt1,
  const vector<double>& pt2)
{
  if (pt1.size() != pt2.size() || pt1.size() != diff.size()) {
    cerr << "Cannot compute vector difference: size mismatch." << endl;
    return;
  }
  for (unsigned i = 0; i < pt1.size(); i++) {
    diff[i] = pt1[i] - pt2[i];
  }
}

VecDbl surfpack::weightedAvg(const VecDbl& first, const VecDbl& second,
  double alpha)
{
  assert(first.size() == second.size());
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  VecDbl result(first.size());
  double beta = 1.0 - alpha;
  for (unsigned i = 0; i < result.size(); i++) {
    result[i] = first[i]*alpha + beta*second[i];
  }
  return result; 
}

// ____________________________________________________________________________
// Functions for common linear algebra tasks 
// ____________________________________________________________________________
void surfpack::linearSystemLeastSquares(MtxDbl& A, 
  vector<double>& x, vector<double> b)
{
  ///\todo Change input b to be pass by value.  It gets overwritten
  /// inside of dgels, and that has caused some problems.  It's probably
  /// better to just copy the whole vector and pass the copy to dgels.
  const int debug = 0;
  if (debug) {
    dbg(debug) << "\nA(" << A.getNRows() <<  "," << A.getNCols() << "):\n";
    dbg(debug) << A.asString() << "\n\n";
    for (unsigned i = 0; i < b.size(); i++) {
      dbg(debug) << b[i] << "\n";
    }
    dbg(debug) << "\n";
  }
  // Rows in A must == size of b
  assert(A.getNRows() == b.size()); 
  // System must be square or over-constrained
  assert(A.getNRows() >= A.getNCols());
  int n_rows = static_cast<int>(A.getNRows());
  int n_cols = static_cast<int>(A.getNCols());
  // Client may supply a "blank" initialized vector for x
  int lwork = n_rows*n_cols * 2;
  vector<double> work(lwork);
  // values must be passed by reference to Fortran, so variables must be 
  // declared for info, nrhs, trans
  int info;
  int nrhs=1;
  char trans = 'N';
  DGELS_F77(&trans,&n_rows,&n_cols,&nrhs,&A(0,0),&n_rows,&b[0],
	    &n_rows,&work[0],&lwork,&info);
  x = b;
  x.resize(n_cols);
  if (debug) {
    for (unsigned i = 0; i < x.size(); i++) {
      dbg(debug) << x[i] << "\n";
    }
  }
}

void surfpack::leastSquaresWithEqualityConstraints(MtxDbl& A, 
  vector<double>& x, vector<double>& c,
  MtxDbl& B, vector<double>& d)
{
  int m = static_cast<int>(A.getNRows());
  int n = static_cast<int>(A.getNCols());
  int p = static_cast<int>(B.getNRows());
  int B_cols = static_cast<int>(B.getNCols());
  assert(B_cols == n);
  assert(p <= n);
  assert(n <= p + m);
  assert(x.size() == static_cast<unsigned>(n));
  int lwork = (m + n + p);
  lwork *= lwork;
  ///\todo Compute optimal blocksize before running dgglse
  vector<double> work(lwork);
  int info = 0;
  //MtxDbl Acopy = A;
  //MtxDbl Bcopy = B;
  //cout << endl;
  //cout << "A" << endl;
  //cout << A.asString() << endl;
  //cout << A.asArrayString() << endl;
  //cout << endl;
  //cout << "B" << endl;
  //cout << B.asString() << endl;
  //cout << B.asArrayString() << endl;
  //cout << "c" << endl;
  //copy(c.begin(),c.end(),ostream_iterator<double>(cout,"\n"));
  //cout << "d" << endl;
  //copy(d.begin(),d.end(),ostream_iterator<double>(cout,"\n"));
  //cout << "x before" << endl;
  //copy(x.begin(),x.end(),ostream_iterator<double>(cout,"\n"));
  DGGLSE_F77(&m,&n,&p,&A(0,0),&m,&B(0,0),&p,&c[0],&d[0],&x[0],&work[0],&lwork,
	     &info);
  //cout << "x after" << endl;
  //copy(x.begin(),x.end(),ostream_iterator<double>(cout,"\n"));
  //vector<double> result;
  //cout << "B*x after" << endl;
  //surfpack::matrixVectorMult(result,Bcopy,x);
  //copy(result.begin(),result.end(),ostream_iterator<double>(cout,"\n"));
  if (info != 0) throw string("Error in dgglse\n");
}

MtxDbl& surfpack::inverse(SurfpackMatrix< double>& matrix)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  /// \todo compute optimal blocksize
  int lwork = n_cols;  // should be optimal blocksize
  vector<int> ipvt(n_rows);
  vector<double> work(lwork);
  int lda = n_cols;
  int info = 0;
  DGETRF_F77(&n_rows,&n_cols,&matrix(0,0),&lda,&ipvt[0],&info);
  DGETRI_F77(&n_rows,&matrix(0,0),&lda,&ipvt[0],&work[0],&lwork,&info);
  return matrix;
}

MtxDbl& surfpack::LUFact(MtxDbl& matrix,
  vector<int>& ipvt)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  // dgetrf may reorder the rows, the mapping between old and new rows
  // is returned in ipvt.
  ipvt.resize(n_rows);
  int lda = n_cols;
  int info = 0;
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRF_F77(&n_rows,&n_cols,&matrix(0,0),&lda,&ipvt[0],&info);
  //std::cout << "Done with dgetrf" << std::endl;
  return matrix;
}

MtxDbl& surfpack::inverseAfterLUFact(MtxDbl& matrix, vector<int>& ipvt)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  int lwork = n_cols;  // should be optimal blocksize
  vector<double> work(lwork);
  int lda = n_rows;
  int info = 0;
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRI_F77(&n_rows,&matrix(0,0),&lda,&ipvt[0],&work[0],&lwork,&info);
  //std::cout << "Done with getri" << std::endl;
  return matrix;
}

VecDbl& surfpack::matrixVectorMult(VecDbl& result,
  MtxDbl& matrix, VecDbl& the_vector, char trans)
{
  assert((trans == 'N' && matrix.getNCols() == the_vector.size()) ||
	 (trans == 'T' && matrix.getNRows() == the_vector.size()));
  int result_size = (trans == 'N') ? matrix.getNRows() : matrix.getNCols();
  result.resize(result_size);
  //char transpose = 'N';
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  int incx = 1;
  int incy = 1;
  double alpha = 1.0;
  double beta = 0.0;
  DGEMV_F77(&trans,&n_rows,&n_cols,&alpha,&matrix(0,0),&n_rows,&the_vector[0],
	    &incx, &beta, &result[0], &incy);
  return result;
}

MtxDbl& surfpack::matrixMatrixMult(MtxDbl& result, MtxDbl& matrixA, 
  MtxDbl& matrixB, char transA, char transB)
{
  assert(matrixA.getNCols(transA) == matrixB.getNRows(transB));
  result.reshape(matrixA.getNRows(transA),matrixB.getNCols(transB));
  int n_rows = static_cast<int>(matrixA.getNRows(transA));
  int n_cols = static_cast<int>(matrixB.getNCols(transB));
  int k = static_cast<int>(matrixA.getNCols(transA));
  int lda = matrixA.getNRows(); // leading dimension irrespective of trans
  int ldb = matrixB.getNRows(); // leading dimension irrespective of trans
  int ldc = result.getNRows(); // leading dimension irrespective of trans
  double alpha = 1.0; // coefficent of matrix multiply
  double beta = 0.0;  // coefficient of result matrix
  DGEMM_F77(&transA,&transB,&n_rows,&n_cols,&k,&alpha,&matrixA(0,0),&lda,&matrixB(0,0),&ldb,&beta,&result(0,0),&ldc);
  return result;
}
  /// matrix-matrix addition
MtxDbl& surfpack::matrixSum(MtxDbl& result, MtxDbl& matrixA, MtxDbl& matrixB)
{
  assert (matrixA.getNRows() == matrixB.getNRows());
  assert (matrixA.getNCols() == matrixB.getNCols());
  result.reshape(matrixA.getNRows(),matrixA.getNCols());
  /// todo check for self in matrix oper= and then create a += func
  /// todo Then this can be simplified/optimized
  for (unsigned i = 0; i < matrixA.getNRows(); i++) {
    for (unsigned j = 0; j < matrixA.getNCols(); j++) {
      result(i,j) = matrixA(i,j) + matrixB(i,j);
    }
  }
  return result;
}

MtxDbl& surfpack::matrixSubtraction(MtxDbl& result, MtxDbl& matrixA, 
  MtxDbl& matrixB)
{
  assert (matrixA.getNRows() == matrixB.getNRows());
  assert (matrixA.getNCols() == matrixB.getNCols());
  result.reshape(matrixA.getNRows(),matrixA.getNCols());
  /// todo check for self in matrix oper= and then create a += func
  /// todo Then this can be simplified/optimized
  for (unsigned i = 0; i < matrixA.getNRows(); i++) {
    for (unsigned j = 0; j < matrixA.getNCols(); j++) {
      result(i,j) = matrixA(i,j) - matrixB(i,j);
    }
  }
  return result;
}

double surfpack::dot_product(const VecDbl& vector_a, 
	     const VecDbl& vector_b)
{
  assert(vector_a.size() == vector_b.size()); 
  int size = static_cast<int>(vector_a.size());
  int inc = 1;
  // ddot will not violate the constness
  return DDOT_F77(&size, const_cast<double*>( &vector_a[0] ), &inc,
		         const_cast<double*>( &vector_b[0] ), &inc);
}

VecDbl& surfpack::vectorShift(VecDbl& the_vector, double shift_value)
{
  for (VecDbl::iterator iter = the_vector.begin();
       iter != the_vector.end(); ++iter) {
    *iter -= shift_value;
  }
  return the_vector;
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

/// Return the value of the test function specified by parameter name at the 
/// point specified by parameter pt
/// \todo Change if-else construct to lookup in STL map of <name, function>
/// pairs.
double surfpack::testFunction(const string name, const vector<double>& pt)
{
  if (name == "rosenbrock") {
    return surfpack::rosenbrock(pt);
  } else if (name == "sphere") {
    return surfpack::sphere(pt);
  } else if (name == "sumofall") {
    return surfpack::sumofall(pt);
  } else if (name == "simplepoly") {
    return surfpack::simplepoly(pt);
  } else if (name == "moderatepoly") {
    return surfpack::moderatepoly(pt);
  } else if (name == "sinewave") {
    return surfpack::sinewave(pt);
  } else if (name == "quasisine") {
    return surfpack::quasisine(pt);
  } else if (name == "xplussinex") {
    return surfpack::xplussinex(pt);
  } else if (name == "noise") {
    return surfpack::noise(pt);
  } else {
    return surfpack::rastrigin(pt);
  }
}

/// Non-trivial polynomial function
double surfpack::moderatepoly(const vector<double>& pt)
{
  double result = -3.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    switch (i % 3) {
      case 0: result -= 2.0*(x-3.0); break;
      case 1: result += 1.0*(x+3.0)*(x+3.0); break;
      case 2: result += 2.0*(x-3.0)*(pt[(i+2)%3]); break;
    }
  }
  return result;
}

/// Tony Giunta's test function: a little wave mixed with a big wave, plus 
/// noise
double surfpack::quasisine(const vector<double>& pt)
{
  double result = 0.0;
  double c = 16.0/15.0;
  double e = 1.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += sin(c*x-e) + sin(c*x-e)*sin(c*x-e) + .02*sin(40.0*(c*x-e));
  }
  return result;
}

/// f(x) = sigma{i=1 to n}(x_i^2 - 10*cos(2*pi*x) + 10).  With side 
/// constraints near zero (e.g. +/-10) along each dimension, the function 
/// appears highly-multimodal.  With larger bounds, the bowl shape becomes
/// more dominant and the waviness is reduced to noise.
double surfpack::rastrigin(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x-10*cos(4.0*acos(0.0)*x)+10.0;
  }
  return result;
}

/// A multi-dimensional extension of the classic Rosenbrock test function
double surfpack::rosenbrock(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size() - 1; i++) {
    double x = pt[i];
    double xp = pt[i+1];
    result += 100.0*(xp-x*x)*(xp-x*x)+(x-1.0)*(x-1.0); 
  }
  return result;
}

/// f(x) = 3 + sigma{i=1 to n}(2*x_i)
double surfpack::simplepoly(const vector<double>& pt)
{
  double result = 3.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += 2.0*x;
  }
  return result;
}

/// Sum of the sine function along each dimension
double surfpack::sinewave(const vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += sin(x);
  }
  return result;
}

/// Sum of squares along each dimension
double surfpack::sphere(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x;
  }
  return result;
}

/// f(x) = sigma{i=1 to i}(x_i)
double surfpack::sumofall(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    result += pt[i];
  }
  return result;
}

/// f(x) = sigma{i=1 to n}(x_i + sin x_i)
double surfpack::xplussinex(const vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x + sin(x);
  }
  return result;
}

/// Random (different queries for the same point will give different results)
double surfpack::noise(const vector<double>& pt)
{
  // original
  //return static_cast<double>(rand());

  // intent of original?
  //return static_cast<double>(surfpack::shared_rng().randInt());

  // guaranteed same range as original -- may not be right
  return static_cast<double>(surfpack::shared_rng().randInt(RAND_MAX));
}

void surfpack::stripQuotes(std::string& str) 
{
    int pos;
    while ( (pos = str.find('\'')) != std::string::npos) {
      str.erase(pos,pos+1);
    }
}
