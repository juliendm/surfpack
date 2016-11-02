/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURF_DATA_H__
#define __SURF_DATA_H__

#include "surfpack_system_headers.h"
#include "SurfPoint.h"

/// Contains a set of SurfPoint objects.  Supports exclusion of some
/// of the SurfPoint objects, so clients may operate on some
/// subset of the data if they so choose (e.g., in a cross-validation 
/// algorithm).  Contains methods for I/O support.  Does not allow duplicate
/// points.
/// \todo Allow the points to be weighted differently, as would be needed
/// in weighted regression.
class SurfData
{

public:

/// Nested exception class. A bad_surf_data exception is thrown whenever a 
/// client attempts to:
/// 1) Add a SurfPoint to the SurfData object that has a different number of
///    dimensions and/or a different number of response values than another
///    SurfPoint already in the set;
/// 2) invoke addResponse(...) when there are not yet any SurfPoints;
/// 3) invoke addResponse(...) with the wrong number of values;
/// 4) invoke addResponse(...) on an object where the logical size of the 
///    data set does not match the phyiscal size (i.e., some of the SurfPoints
///    have been marked for exclusion);
/// 5) write a SurfData object that contains no data to a stream;
/// 6) create a Surface with a SurfData object that does not have enough points. 
class bad_surf_data : public std::runtime_error
{
public:
  bad_surf_data(const std::string& msg = "") : std::runtime_error(msg) {}
};


// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
  /// Vector of points will be copied and checked for duplicates
  SurfData(const std::vector<SurfPoint>& points_);

  /// Read a set of SurfPoints from a file; binary or text determined
  /// based on extension
  SurfData(const std::string filename);

  /// Read a set of SurfPoints from a std::istream, binary or text as
  /// specified via flag
  SurfData(std::istream& is, bool binary = false);

  /// Read a set of SurfPoints in text format from a std::istream.
  /// The stream does not
  /// contain the normal header information (#points, #vars, #responses).
  /// The #vars and #responses are explicitly specified in the constructor;
  /// The stream reader processes data until eof, assuming one point per line.

  // TODO: include grad/hess sizing?

  SurfData(const std::string filename, unsigned n_vars, unsigned n_responses, 
    unsigned n_cols_to_skip);

  /// Make a deep copy of the object 
  SurfData(const SurfData& other); 
  
  /// STL data members' resources automatically deallocated 
  ~SurfData();

  /// Copy only the points which have not been marked for exclusion
  SurfData copyActive();
  
  /// First SurfPoint added will determine the dimensions of the data set 
  SurfData();

private:
  /// Data member initialization that is common to all constructors
  void init();
  
  /// Call delete on the SurfPoint* in the data set.
  void cleanup();

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:
  /// Makes deep copy 
  SurfData& operator=(const SurfData& other);

  /// Makes deep comparison
  bool operator==(const SurfData& other) const;

  /// Makes deep comparison
  bool operator!=(const SurfData& other) const;

  /// Return a const reference to active SurfPoint at given index
  const SurfPoint& operator[](unsigned index) const;

  /// Return the x-value for point pt along dimension dim
  double operator()(unsigned pt, unsigned dim) const;

  /// Return the vector of predictor (x) vars for point index 
  const std::vector<double>& operator()(unsigned pt) const;

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

  /// Return the number of active SurfPoints in the data set 
  unsigned size() const;

  /// True if there are no points
  bool empty() const;
  
  /// Return the dimensionality of the SurfPoints 
  unsigned xSize() const;

  /// Return the number of response functions in the data set
  unsigned fSize() const;

  /// Return the set of excluded points (the indices)
  const std::set<unsigned>& getExcludedPoints() const ; 

  /// Get the default response f value of the (index)th point
  double getResponse(unsigned index) const;

  /// Get the vector of gradients for the (index)th point
  //  const std::vector<double>& getGradient(unsigned index) const;

  /// Get the full symmetric matrix of Hessians for the (index)th point
  //  const SurfpackMatrix<double>& getHessian(unsigned index) const;

  /// Get the default response f values for all of the points as a vector
  std::vector< double > getResponses() const;
  
  /// Get the predictor (index-th x-value) for all the active points as a vector
  std::vector< double > getPredictor(unsigned index) const;

  /// Get the constraint point
  const SurfPoint& getConstraintPoint() const;

  /// Get the number of constraint data included in the provided point
  unsigned numConstraints() const;

  /// Return defaultIndex
  unsigned getDefaultIndex() const;

  /// Retrieve the label for one of the predictor variables
  const std::string& getXLabel(unsigned index) const;

  /// Retrieve the label for one of the response variables
  const std::string& getFLabel(unsigned index) const;

  /// Retrieve the index and variable type (predictor/response) for a given
  /// name.  Return false if not found
  bool varIndex(const std::string& name, unsigned& index, bool& isResponse) const;

  

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Specify which response value getResponse will return. When a Surface 
  /// object that is associated with the SurfData object operates on the data,
  /// it sets this value so that the response value lookup function will return
  /// the value for the response variable that that particular Surface object
  /// is interested in.  
  void setDefaultIndex(unsigned index) const; 
  
  /// Set the response value of the (index)th point that corresponds to this
  /// surface
  void setResponse(unsigned index, double value);

  /// Add a point to the data set. The passed point will be copied.
  void addPoint(const SurfPoint& sp);

  /// Add a new response variable to each point, possibly with label. 
  /// Return the index of the new variable.
  unsigned addResponse(const std::vector<double>& newValues, 
		       std::string label = ""); 

  /// Set the constraint (anchor) point, copying the data (only single
  /// constraint supported)
  void setConstraintPoint(const SurfPoint& sp);
  
  /// Specify which points should be skipped.  This can be used when only a 
  /// subset of the SurfPoints should be used for some computation.
  void setExcludedPoints(const std::set<unsigned>& excluded_points);

  /// For use with copy constructor and assignment operator-- creates a list of
  /// pointers to the points in the data set which is used to check for 
  /// duplication when other points are added in the future
  void buildOrderedPoints();

  /// Set the labels for the predictor variables
  void setXLabels(const std::vector<std::string>& labels);

  /// Set the labels for the response variables
  void setFLabels(const std::vector<std::string>& labels);

  /// Set the label for a single response variable
  void setFLabel(unsigned index, const std::string& label);

private:
  /// Maps all indices to themselves in the mapping data member, based
  /// on number of points present when called.
  void defaultMapping();

  /// Set x vars labels to 'x0' 'x1', etc.; resp. vars to 'f0' 'f1', etc.
  void defaultLabels();

public:
   
// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write the active SurfPoints to a file
  /// binary extension .bspd: includes header not labels
  /// text extension    .spd: includes header and label info
  /// text extension    .dat: no header or labels
  void write(const std::string& filename) const;

  /// Read a set of SurfPoints from a file
  void read(const std::string& filename);
  
  /// Write the data in binary format, with header, but no labels
  void writeBinary(std::ostream& os) const;

  /// Write the data in text format, with optional header, optional
  /// labels, followed by points.
  void writeText(std::ostream& os, bool write_header = true,
		 bool write_labels = true) const; 

  /// Read the data in binary format
  void readBinary(std::istream& is); 

  /// Read the data in text format
  void readText(std::istream& is, bool read_header = true, 
    unsigned skip_columns = 0); 

private:

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

  /// Dimensionality of the input (x) space from which the SurfPoints are drawn
  unsigned xsize;

  /// Number of response/output variables in the data set 
  unsigned fsize;

  /// Number of responses with gradient data (must be 0 or fsize)
  unsigned gradsize;

  /// Number of responses with Hessian data (must be 0 or fsize)
  unsigned hesssize;

  /// The set of points in this data set
  std::vector<SurfPoint*> points; 

  /// The indices of points that are to be excluded in computation. This can
  /// be used in a cross-validation scheme to systematically ignore parts of
  /// data set at different times.  
  std::set<unsigned> excludedPoints;

  /// For mapping the indices in points to the indices returned by operator[].
  /// Normally, mapping[i] is equal to i, but the set of excludedPoints is not
  /// empty, this will not be the case.
  std::vector<unsigned> mapping;

  /// The index of the response variable that will be returned by F
  mutable unsigned defaultIndex;

  /// Constraint (anchor) point which the model should match exactly
  /// point). Some models require separate treatment, e.g.,
  /// LinearRegressionModel, though KrigingModel would allow
  /// integrated in the points array.  For other surrogates, these are
  /// typically just integrated into the points.
  SurfPoint constraintPoint;

  /// Labels for the predictor/input (x) variables
  std::vector< std::string > xLabels;
 
  /// Labels for the responses/output variables
  std::vector< std::string > fLabels;

public:
  typedef std::set<SurfPoint*,SurfPoint::SurfPointPtrLessThan> SurfPointSet;

private:
  /// Stores the same set of SurfPoint* that points does, but because it is a 
  /// set, membership tests can be done in O(log n) time.  When combined with
  /// the SurfPointPtrLessThan functor object, it allows a SurfData object to
  /// check all SurfPoints in the data set against all others for duplication.
  /// This can be done in O(n log n) time instead of O(n^2).
  SurfPointSet orderedPoints;

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  /// Returns true if file has .bspd extension, false if it has .spd extension. 
  /// Otherwise, an exception is thrown.
  bool hasBinaryFileExtension(const std::string& filename) const;

  /// If the line contains single-quoted string, parse them out as labels
  /// and return true; otherwise, return false
  bool readLabelsIfPresent(std::string single_line);

  /// Read the #points, #vars, #responses
  unsigned readHeaderInfo(std::istream& is);

  /// return the sample points (x data) as a matrix samples[i][j],
  /// where i indexes sample and j indexes dimension of x
  static VecVecDbl asVecVecDbl(const SurfData& data);

      
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class SurfData data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif
  

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

protected:
  // Throw an exception if there are any mismatches in the number of
  // dimensions or number of response values among points in the data set  
  void sanityCheck() const;

  /// Make sure an index falls within acceptable boundaries
  void checkRangeNumPoints(const std::string& header, unsigned index) const;

  /// Make sure an index falls within acceptable boundaries
  void checkRangeNumResponses(const std::string& header, unsigned index) const;

#ifdef __TESTING_MODE__ 
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

/// Print the SurfData to a stream 
std::ostream& operator<<(std::ostream& os, const SurfData& data);


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SurfData::serialize(Archive & archive, 
			 const unsigned int version)
{  
  archive & xsize;
  archive & fsize;
  archive & gradsize;
  archive & hesssize;
  archive & points;
  archive & excludedPoints;
  archive & mapping;
  archive & defaultIndex;
  archive & constraintPoint;
  archive & xLabels;
  archive & fLabels;
  archive & orderedPoints;
}
#endif



#endif
