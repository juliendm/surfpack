/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "SurfPoint.h"

using std::cerr;
using std::endl;
using std::ios;
using std::istream;
using std::istringstream;
using std::ostream;
using std::ostringstream;
using std::range_error;
using std::setw;
using std::string;
using std::vector;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(SurfPoint)
#endif


// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Initialize without any response values
SurfPoint::SurfPoint(const vector<double>& x) : x(x)
{
  init();
}

/// Initialize point with one response value
SurfPoint::SurfPoint(const vector<double>& x, double f0) : x(x), f(1)
{
  f[0] = f0;
  init();
}

/// Initialize with one response value and corresponding gradient
SurfPoint::SurfPoint(const std::vector<double>& x, double f0, 
		     const std::vector<double>& gradient0) 
  : x(x), f(1), fGradients(1)
{
  f[0] = f0;
  fGradients[0] = gradient0;
  init();
}

/// Initialize with one response value and corresponding gradient and Hessian
SurfPoint::SurfPoint(const std::vector<double>& x, double f0, 
		     const std::vector<double>& gradient0,
		     const SurfpackMatrix<double>& hessian0)
  : x(x), f(1), fGradients(1), fHessians(1)
{
  f[0] = f0;
  fGradients[0] = gradient0;
  fHessians[0] = hessian0;
  init();
}

/// Initialize with zero or more response values
SurfPoint::SurfPoint(const vector<double>& x, const vector<double>& f)
  : x(x), f(f)
{
  init();
}

/// Read point from istream in binary format
SurfPoint::SurfPoint(istream& is, unsigned xsize, unsigned fsize, 
		     unsigned grad_size, unsigned hess_size) 
  : x(xsize), f(fsize), fGradients(grad_size), fHessians(hess_size)
{
  for(unsigned i = 0; i < fsize; ++i) {
    fGradients.resize(xsize);
    fHessians.resize(xsize, xsize);
  }
  readBinary(is);
  init();
}

/// Read point from string in text format
SurfPoint::SurfPoint(const string& single_line, unsigned xsize, unsigned fsize, 
		     unsigned grad_size, unsigned hess_size,
		     unsigned skip_columns) 
  : x(xsize), f(fsize)
{
  readText(single_line, skip_columns);
  init();
}

/// Copy constructor performs a deep copy
SurfPoint::SurfPoint(const SurfPoint& sp) 
  : x(sp.x), f(sp.f), fGradients(sp.fGradients), fHessians(sp.fHessians)
{
  init();
}


/// Default constructor creates a one dimensional point at the origin 
/// with no function value (empty)
SurfPoint::SurfPoint() : x(1), f(0)
{
  x[0] = 0;
  init();
} 

/// Initialization used by all regular constructors.  Ensures that point has
/// at least one dimension.
void SurfPoint::init()
{
  if (x.empty()) {
    throw SurfPoint::null_point();
  }
  // if any gradient data, must have for all functions
  if (!fGradients.empty() && f.size() != fGradients.size()) {
    throw SurfPoint::bad_gradient_size();
  }
  // if any Hessian data, must have for all functions
  if (!fHessians.empty() && f.size() != fHessians.size()) {
    throw SurfPoint::bad_gradient_size();
  }
}

SurfPoint::~SurfPoint() 
{

}

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// Assign 'other' to 'this' unless they are already equal 
SurfPoint& SurfPoint::operator=(const SurfPoint& other) 
{
  if (*this != other) {
    x = other.x;
    f = other.f;
    fGradients = other.fGradients;
    fHessians = other.fHessians;
  }
  return (*this);
}

bool doubles_match(double x, double y)
{
    if (fabs(x) < 1e-10) {
      if (fabs(y) > 1e-10) return false;
    } else {
      if (fabs(x-y)/fabs(x) > 1e-10) return false;
    }
    return true;
}
/// Tests for deep equality
bool SurfPoint::operator==(const SurfPoint& other) const
{
  //return x == other.x && f == other.f;
  for (unsigned i = 0; i < x.size(); i++) {
    if (!doubles_match(x[i],other.x[i])) return false;
  }
  for (unsigned i = 0; i < f.size(); i++) {
    if (!doubles_match(f[i],other.f[i])) return false;
  }
  // gradients
  for (unsigned fi = 0; fi < fGradients.size(); ++fi) {
    for (unsigned xj = 0; xj < xSize(); ++xj) {
      if (!doubles_match(fGradients[fi][xj], other.fGradients[fi][xj])) 
	return false;
    }
  }
  // Hessians
  for (unsigned fi = 0; fi < fHessians.size(); ++fi) {
    for (unsigned xj = 0; xj < xSize(); ++xj) {
      for (unsigned xk = 0; xk < xSize(); ++xk) {
	if (!doubles_match(fHessians[fi](xj,xk), other.fHessians[fi](xj,xk)))
	  return false;
      }
    }
  }
  return true;
}

/// Tests for deep inequality
bool SurfPoint::operator!=(const SurfPoint& other) const
{
  return !(*this == other);
} 

/// Return the value along the (xindex)th dimension;
double SurfPoint::operator[](unsigned xindex) const
{
  if (xindex >= x.size()) throw string("Out of range in SurfPoint");
  return x[xindex];
}
  
/// Function object for use with pairs of SurfPoint objects (particularly in
/// a SurfData object).  SurfPoint s1 is "less than" s2 if it has fewer
/// dimensions.  SurfPoint s1 is also less than s2 if, for some dimension i,
/// s1[i] < s2[i] AND s1[j] == s2[j] for all j, 0 <= j < i.  Note that the
/// SurfPoint's response values have no bearing on the results for this 
/// comparison.  Since the response values DO affect the results of 
/// SurfPoint::operator==, it is NOT necessarily the case that s1 == s2 and
/// (!SurfPointPtrLessThan(&s1,&s2) && !SurfPointPtrLessThan(&s2,&s1)) will
/// return the same boolean value.
bool SurfPoint::SurfPointPtrLessThan::
  operator()(const SurfPoint* sp1, const SurfPoint* sp2) const
{
  if (sp1->X().size() < sp2->X().size()) {
    return true;
  } else if (sp1->X().size() > sp2->X().size()) {
    return false;
  } else {
    for (unsigned i = 0; i < sp1->X().size(); i++) {
      if (sp1->X()[i] < sp2->X()[i]) {
        return true;
      } else if (sp1->X()[i] > sp2->X()[i]) {
        return false;
      }
    }
    // If execution makes it this far, the points have the
    // same dimensionality and are equivalent along each dimension.
    return false;
  }
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return dimensionality of data point
unsigned SurfPoint::xSize() const
{ 
  return x.size(); 
}

/// Return number of response variables
unsigned SurfPoint::fSize() const
{ 
  return f.size(); 
}

/// Return number of response variables having gradients
unsigned SurfPoint::fGradientsSize() const
{ 
  return fGradients.size(); 
}

/// Return number of response variables
unsigned SurfPoint::fHessiansSize() const
{ 
  return fHessians.size(); 
}

/// Return point in the domain 
const vector<double>& SurfPoint::X() const
{ 
  return x; 
}

/// Return response value at responseIndex
double SurfPoint::F(unsigned responseIndex) const
{ 
  string header(
    "Error in query SurfPoint::F. Invalid responseIndex."
  );
  // Throw an exception if the responseIndex is out of range.
  checkRange(header, responseIndex);
  return f[responseIndex]; 
}

/// Return gradient vector for response value at responseIndex
const std::vector<double>& SurfPoint::fGradient(unsigned responseIndex) const
{
  string header(
    "Error in query SurfPoint::fGradient. Invalid responseIndex."
  );
  // Throw an exception if the responseIndex is out of range.
  checkRange(header, responseIndex);
  return fGradients[responseIndex];
}

/// Return response value at responseIndex
const SurfpackMatrix<double>& SurfPoint::fHessian(unsigned responseIndex) const
{
  string header(
    "Error in query SurfPoint::fHessian. Invalid responseIndex."
  );
  // Throw an exception if the responseIndex is out of range.
  checkRange(header, responseIndex);
  return fHessians[responseIndex];
}


// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Append a new response variable
unsigned SurfPoint::addResponse(double val)
{
  f.push_back(val);
  return f.size()-1;
}

/// Set an existing response variable to a new value
void SurfPoint::F(unsigned responseIndex, double responseValue)
{ 
  string header(
    "Error in command SurfPoint::F. Invalid responseIndex. No update made."
  );
  // Throw an exception if the responseIndex is out of range.
  checkRange(header, responseIndex);
  f[responseIndex] = responseValue; 
}

/// Change the value of one of the dimensions of the point
void SurfPoint::setX(unsigned index, double value)
{
  if (index >= x.size()) {
    x.resize(index+1);
  }
  x[index] = value;
}

/// Change the dimensionality of the point
void SurfPoint::resize(unsigned new_size)
{
  x.resize(new_size);
}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write location (x) and responses (f) of this point to stream in
/// binary format
void SurfPoint::writeBinary(ostream& os) const
{
  // write x
  for (unsigned i = 0; i < x.size(); i++) {
    os.write(reinterpret_cast<const char*>(&x[i]),sizeof(x[i])) ;
  }
  // write f
  for (unsigned i = 0; i < f.size(); i++) {
    os.write(reinterpret_cast<const char*>(&f[i]),sizeof(f[i]));
  }
  // write gradient
  for (unsigned i = 0; i < fGradients.size(); ++i) {
    for (unsigned j = 0; j < x.size(); ++j) {
      os.write(reinterpret_cast<const char*>(&fGradients[i][j]),
	       sizeof(fGradients[i][j]));
    }
  }
  // write Hessian
  for (unsigned i = 0; i < fHessians.size(); ++i) {
    for (unsigned j = 0; j < x.size(); ++j) {
      for (unsigned k = 0; k < x.size(); ++k) {
	os.write(reinterpret_cast<const char*>(&fHessians[i](j, k)),
		 sizeof(fHessians[i](j, k)));
      }
    }
  }
}

/// Write location and responses of this point to stream in text format 
void SurfPoint::writeText(ostream& os) const
{
  // ios_base::flags returns ios::fmtflags object, but OSF compiler doesn't 
  // like that.
  // Save the stream flags.  The output precision may be modified, but it 
  // will be restored to its old value before the method exits. 
  ios::fmtflags old_flags = os.flags();
  unsigned old_precision = os.precision(surfpack::output_precision);
  os.setf(ios::scientific);
  // write x
  for (unsigned i = 0; i < x.size(); ++i) {
    os  << setw(surfpack::field_width) << x[i] ;
  }
  // write f
  for (unsigned i = 0; i < f.size(); ++i) {
    os << setw(surfpack::field_width) << f[i];
  }
  // write gradient
  for (unsigned i = 0; i < fGradients.size(); ++i) {
    for (unsigned j = 0; j < x.size(); ++j) {
      os << setw(surfpack::field_width) << fGradients[i][j];
    }
  }
  // write Hessian
  for (unsigned i = 0; i < fHessians.size(); ++i) {
    for (unsigned j = 0; j < x.size(); ++j) {
      for (unsigned k = 0; k < x.size(); ++k) {
	os << setw(surfpack::field_width) << fHessians[i](j, k);
      }
    }
  }
  os << endl;
  // Restore output flags to what they were before this method was called.
  os.flags(old_flags);
  os.precision(old_precision);
}
 
void SurfPoint::readBinary(istream& is)
{
  unsigned xValsRead = 0;
  unsigned fValsRead = 0;
  unsigned fGradRead = 0;
  unsigned fHessRead = 0;
  try {
    // read the point in binary format
    // For each read, throw an exception via checkForEOF if there are
    // fewer values on the line than expected.
    // read x
    for (xValsRead = 0; xValsRead < x.size(); xValsRead++) {
       surfpack::checkForEOF(is);
       is.read(reinterpret_cast<char*>(&x[xValsRead]),sizeof(x[xValsRead]));
    }
    // read f
    for (fValsRead = 0; fValsRead < f.size(); fValsRead++) {
       surfpack::checkForEOF(is);
       is.read(reinterpret_cast<char*>(&f[fValsRead]),sizeof(f[fValsRead]));
    }
    // read gradient
    for (fGradRead = 0; fGradRead < fGradients.size(); ++fGradRead) {
      for (unsigned j = 0; j < x.size(); ++j) {
	surfpack::checkForEOF(is);
	is.read(reinterpret_cast<char*>(&fGradients[fGradRead][j]),
		sizeof(fGradients[fGradRead][j]));
      }
    }
    // read Hessian
    for (fHessRead = 0; fHessRead < fHessians.size(); ++fHessRead) {
      for (unsigned j = 0; j < x.size(); ++j) {
	for (unsigned k = 0; k < x.size(); ++k) {
	  surfpack::checkForEOF(is);
	  is.read(reinterpret_cast<char*>(&fHessians[fHessRead](j,k)),
		  sizeof(fHessians[fHessRead](j,k)));
	}
      }
    }
  } catch (surfpack::io_exception&) {
    cerr << "Bad binary read. "
	 << "\nExpected on this line: " << x.size() << " domain value(s) "
         << "and " << f.size() << " response value(s);" << endl;
    if (fGradients.size() > 0 || fHessians.size() > 0)
      cerr << fGradients.size() << " gradient(s) and "
	   << fHessians.size() << " Hessians(s);" << endl;
    cerr << "Found: " << xValsRead << " domain value(s), " 
	 << fValsRead << " response value(s),\n" 
	 << fGradRead << " gradient(s), and "
	 << fHessRead << " Hessian(s)." << endl;
    throw;
  } catch (...) {
    cerr << "Exception rethrown in SurfPoint::readBinary(istream& is)" 
         << endl;
    throw;
  }
}

void SurfPoint::readText(const string& single_line, unsigned skip_columns) 
{
  unsigned xValsRead = 0;
  unsigned fValsRead = 0;
  unsigned fGradRead = 0;
  unsigned fHessRead = 0;
  string dummy;
  try {
    // read the point as text
    istringstream streamline(single_line);
    for (unsigned i = 0; i < skip_columns; i++) streamline >> dummy;
    // For each read, throw an exception via checkForEOF if there are
    // fewer values on the line than expected.
    // read x
    for (xValsRead = 0; xValsRead < x.size(); ++xValsRead) {
       surfpack::checkForEOF(streamline);
       streamline >> x[xValsRead];
    }
    // read f
    for (fValsRead = 0; fValsRead < f.size(); ++fValsRead) {
       surfpack::checkForEOF(streamline);
       streamline >> f[fValsRead];
    }
    // read gradient
    for (fGradRead = 0; fGradRead < fGradients.size(); ++fGradRead) {
      for (unsigned j = 0; j < x.size(); ++j) {
	surfpack::checkForEOF(streamline);
	streamline >> fGradients[fGradRead][j];
      }
    }
    // read Hessian
    for (fHessRead = 0; fHessRead < fHessians.size(); ++fHessRead) {
      for (unsigned j = 0; j < x.size(); ++j) {
	for (unsigned k = 0; k < x.size(); ++k) {
	  surfpack::checkForEOF(streamline);
	  streamline >> fHessians[fGradRead](j,k);
	}
      }
    }
  } catch (surfpack::io_exception&) {
    cerr << "Bad SurfPoint: " << single_line 
	 << "\nExpected on this line: " << x.size() << " domain value(s) "
         << "and " << f.size() << " response value(s);" << endl;
    if (fGradients.size() > 0 || fHessians.size() > 0)
      cerr << fGradients.size() << " gradient(s) and "
	   << fHessians.size() << " Hessians(s);" << endl;
    cerr << "Found: " << xValsRead << " domain value(s), " 
	 << fValsRead << " response value(s),\n" 
	 << fGradRead << " gradient(s), and "
	 << fHessRead << " Hessian(s)." << endl;
    throw;
  } catch (...) {
    cerr << "Exception caught and rethrown in SurfPoint::readText(istream& is)"
         << endl;
    throw;
  }
}

/// Write point to an output stream in text format
ostream& operator<<(ostream& os, const SurfPoint& sp) 
{
  sp.writeText(os);
  return os;
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

/// Provides range checking on the response values.  Throws an exception if an
/// index is requested that does not exist.
void SurfPoint::checkRange(const string& header, unsigned index) const
{
  if (index >= f.size()) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (f.empty()) {
      errormsg << "There are no response values associated with this point"
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << f.size() - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}
