/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "SurfData.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ios;
using std::list;
using std::ofstream;
using std::ostream;
using std::ostream_iterator;
using std::ostringstream;
using std::range_error;
using std::set;
using std::setw;
using std::string;
using std::vector;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(SurfData)
#endif


// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Vector of points will be copied and checked for duplicates
SurfData::SurfData(const vector<SurfPoint>& points_) 
{
  if (points_.empty()) {
    this->xsize = 0;
    this->fsize = 0;
    this->gradsize = 0;
    this->hesssize = 0;
  } else {
    this->xsize = points_[0].xSize();
    this->fsize = points_[0].fSize();
    this->gradsize = points_[0].fGradientsSize();
    this->hesssize = points_[0].fHessiansSize();
    defaultLabels();
    for (unsigned i = 0; i < points_.size(); i++) {
      this->addPoint(points_[i]);
    }
  }
  init();
  // Check to make sure data points all have the same number of dimensions
  // and response values.  An exception will be thrown otherwise.
  sanityCheck();
}

/// Read a set of SurfPoints from a file
SurfData::SurfData(const string filename):
  xsize(0), fsize(0), gradsize(0), hesssize(0) 
{
  init();
  read(filename);
}
  
/// Read a set of SurfPoints from a istream.  The stream does not
/// contain the normal header information (#points, #vars, #responses).
/// The #vars and #responses are explicitly specified in the constructor;
/// The stream reader processes data until eof, assuming one point per line.
SurfData::SurfData(const string filename, unsigned n_vars, 
		   unsigned n_responses, unsigned n_cols_to_skip):
  xsize(n_vars), fsize(n_responses), gradsize(0), hesssize(0)
{
  init();

  if (!surfpack::hasExtension(filename,".dat") 
    && !surfpack::hasExtension(filename,".spd")) {
    cerr << "Bad filename: " << filename << endl;
    throw surfpack::io_exception(
      "Expected .dat extension for filename"
    );
  }
  ifstream infile(filename.c_str(),ios::in);
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  }
  bool read_header = false;
  readText(infile, read_header, n_cols_to_skip);
}

/// Read a set of SurfPoints from a istream
SurfData::SurfData(istream& is, bool binary):
  xsize(0), fsize(0), gradsize(0), hesssize(0)
{
  init();
  if (binary) {
    readBinary(is);
  } else {
    readText(is);
  }
}

/// Makes a deep copy of the object 
SurfData::SurfData(const SurfData& other):
  xsize(other.xsize), fsize(other.fsize),
  gradsize(other.gradsize), hesssize(other.hesssize),
  excludedPoints(other.excludedPoints), defaultIndex(other.defaultIndex),
  xLabels(other.xLabels), fLabels(other.fLabels)
{
  for (unsigned i = 0; i < other.points.size(); i++) {
    this->addPoint(*other.points[i]);
  }
  mapping = other.mapping;
  buildOrderedPoints();
}

/// First SurfPoint added will determine the dimensions of the data set 
SurfData::SurfData(): xsize(0), fsize(0), gradsize(0), hesssize(0)
{
    init();
}

/// STL data members' resources automatically deallocated 
SurfData::~SurfData() 
{
  cleanup();
}

/// Data member initialization that is common to all constructors
void SurfData::init()
{
  defaultIndex = 0;
  defaultMapping();
}

/// Copy only the points which have not been marked for exclusion
SurfData SurfData::copyActive()
{
  // This is not a very efficient method.  It is not expected that it will
  // be called often.  The active points are copied into the local
  // activePoints vector.  Then, a new SurfData object is created, which 
  // requires that all the points be recopied.  Then, at the end of the method
  // the new SurfData object is returned by value, which means all of the data
  // are most likely copied again.  If this is a bottleneck, it would not be
  // terribly difficult to eliminate at least two of the three copies.
  vector<SurfPoint> activePoints;
  for (unsigned i = 0; i < mapping.size(); i++) {
    activePoints.push_back(*points[mapping[i]]);
  }
  SurfData newSD(activePoints);
  if (!activePoints.empty()) {
    newSD.setDefaultIndex(defaultIndex);
  }
  return newSD;
}
  
/// Call delete on the SurfPoint* in the data set.
void SurfData::cleanup()
{
  mapping.clear();
  orderedPoints.clear();
  for (unsigned j = 0; j < points.size(); j++) {
    delete points[j];
    points[j] =0;
  }
  points.clear();
  excludedPoints.clear();
}

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// Makes deep copy of other
SurfData& SurfData::operator=(const SurfData& other)
{
  if (*this != other) {
    xLabels = other.xLabels;
    fLabels = other.fLabels;
    cleanup();
    this->xsize = other.xsize;
    this->fsize = other.fsize;
    this->gradsize = other.gradsize;
    this->hesssize = other.hesssize;
    for (unsigned i = 0; i < other.points.size(); i++) {
      this->addPoint(*other.points[i]);
    }
    this->excludedPoints = other.excludedPoints;
    this->mapping = other.mapping;
    this->defaultIndex = other.defaultIndex;
  }
  buildOrderedPoints();
  return (*this);
}

/// Makes deep comparison
bool SurfData::operator==(const SurfData& other) const
{
  if (this->xsize == other.xsize && 
      this->fsize == other.fsize &&
      this->gradsize == other.gradsize &&
      this->hesssize == other.hesssize &&
      this->size() == other.size()) { 
    for (unsigned i = 0; i < points.size(); i++) {
      if (*this->points[i] != *other.points[i]) {
        return false;
      }
    }
    return true;
  } else {
    return false;
  }
}
      
/// Makes deep comparison
bool SurfData::operator!=(const SurfData& other) const
{
  return !(*this == other);
}

/// Return a const reference to SurfPoint at given index
const SurfPoint& SurfData::operator[](unsigned index) const
{
  static string header("Indexing error in SurfData::operator[] const.");
  checkRangeNumPoints(header, index);
  return *points[mapping[index]];
}

/// Return the x-value for point pt along dimension dim
double SurfData::operator()(unsigned pt, unsigned dim) const
{
  assert(pt < size());
  assert(dim < xSize());
  return points[mapping[pt]]->X()[dim];
}

/// Return the vector of predictor vars for point index 
const vector<double>& SurfData::operator()(unsigned pt) const
{
  if (pt >= size()) {
    cout << "Assertion failure.  Pt: " << pt << " size: " << size() << endl;
  }
  assert(pt < size());
  return points[mapping[pt]]->X();
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return the number of SurfPoints in the data set 
unsigned SurfData::size() const 
{ 
  return mapping.size(); 
}

/// True if there are no points
bool SurfData::empty() const
{
  return mapping.empty();
}

/// Return the dimensionality of the SurfPoints 
unsigned SurfData::xSize() const 
{ 
  return xsize; 
}

/// Return the number of response variables in the data set
unsigned SurfData::fSize() const 
{ 
  return fsize; 
}

/// Return the set of excluded points (the indices)
const set<unsigned>& SurfData::getExcludedPoints() const 
{
  return excludedPoints;
}

/// Get the response value of the (index)th point that corresponds to this
/// surface
double SurfData::getResponse(unsigned index) const
{
  static string header("Indexing error in SurfData::getResponse.");
  checkRangeNumPoints(header, index);
  return points[mapping[index]]->F(defaultIndex);
}

// const std::vector<double>& SurfData::getGradient(unsigned index) const
// {
//   static string header("Indexing error in SurfData::getResponse.");
//   checkRangeNumPoints(header, index);
//   return points[mapping[index]]->fGradient(defaultIndex);
// }

// const SurfpackMatrix<double>& SurfData::getHessian(unsigned index) const
// {
//   static string header("Indexing error in SurfData::getResponse.");
//   checkRangeNumPoints(header, index);
//   return points[mapping[index]]->fHessian(defaultIndex);
// }

/// Get the default response for all of the points as a vector
std::vector< double > SurfData::getResponses() const
{
  vector< double > result(mapping.size());
  for (unsigned i = 0; i < mapping.size(); i++) {
    result[i] = points[mapping[i]]->F(defaultIndex);
  }
  return result;
}

/// Get the predictor (index-th x-value) for all the active points as a vector
std::vector< double > SurfData::getPredictor(unsigned index) const
{
  assert(index < xSize());
  vector< double > result(mapping.size());
  for (unsigned i = 0; i < mapping.size(); i++) {
    result[i] = (*this)(i,index);
  }
  return result;
}

const SurfPoint& SurfData::getConstraintPoint() const
{
  return constraintPoint;
}

/** Calculate the number of constraint data present in functions,
    gradients, and Hessians in the constraintPoint. */
unsigned SurfData::numConstraints() const
{
  unsigned num_constraints = 0;
  if (constraintPoint.fSize() > 0)
    num_constraints += 1; // value at a particular point
  if (constraintPoint.fGradientsSize() > 0)
    num_constraints += xsize; // gradient at a point
  if (constraintPoint.fHessiansSize() > 0) 
    num_constraints += (xsize*xsize+xsize)/2; // hessian at a point
  return num_constraints;
}


/// Get default index
unsigned SurfData::getDefaultIndex() const
{
  return defaultIndex;
}

/// Retrieve the label for one of the predictor variables
const string& SurfData::getXLabel(unsigned index) const
{
  return xLabels[index];
}

/// Retrieve the label for one of the response variables
const string& SurfData::getFLabel(unsigned index) const
{
  return fLabels[index];
}

/// Retrieve the index and variable type (predictor/response) for a given
/// name.  Return false if not found
bool SurfData::varIndex(const string& name, unsigned& index, 
  bool& isResponse) const
{
  //\todo check for apostrophes rather than just assuming they're there
  // Strip off the apostrophes at the beginning and end
  string unquoted_name = name;
  if (name.find('\'') != string::npos) {
    unquoted_name = string(name,1,name.size()-2);
  }
  vector< string>::const_iterator iter = 
    find(xLabels.begin(),xLabels.end(),unquoted_name);
  if (iter != xLabels.end()) {
    index = iter - xLabels.begin();
    isResponse = false;
    return true;
  } else {
    iter = find(fLabels.begin(),fLabels.end(),unquoted_name);
    if (iter != fLabels.end()) {
      index = iter - fLabels.begin();
      isResponse = true;
      return true;
    } else {
      cout << "Name sought: " << unquoted_name << endl;
      cout << "Predictors: " << endl;
      copy(xLabels.begin(),xLabels.end(),ostream_iterator<string>(cout,"\n"));
      cout << "Responses: " << endl;
      copy(fLabels.begin(),fLabels.end(),ostream_iterator<string>(cout,"\n"));
      return false;
    }
  }
}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Specify which response value getResponse will return. When a Surface 
/// object that is associated with the SurfData object operates on the data,
/// it sets this value so that the response value lookup function will return
/// the value for the response variable that that particular Surface object
/// is interested in.  
void SurfData::setDefaultIndex(unsigned index) const
{
  static string header("Indexing error in SurfData::setDefaultIndex.");
  checkRangeNumResponses(header, index);
  defaultIndex = index;
}
  
/// Set the response value of the (index)th point that corresponds to this
/// surface
void SurfData::setResponse(unsigned index, double value)
{
  static string header("Indexing error in SurfData::setResponse.");
  checkRangeNumPoints(header, index);
  points[mapping[index]]->F(defaultIndex, value);
}
  
/// Add a point to the data set. The parameter point will be copied.
void SurfData::addPoint(const SurfPoint& sp) 
{
  if (points.empty()) {
    xsize = sp.xSize();
    fsize = sp.fSize();
    gradsize = sp.fGradientsSize();
    hesssize = sp.fHessiansSize();
    if (xLabels.empty()) {
      defaultLabels();
    }
  } else {
    if (sp.xSize() != xsize || sp.fSize() != fsize ||
	sp.fGradientsSize() != gradsize || sp.fHessiansSize() != hesssize) {
      ostringstream errormsg;
      errormsg << "Error in SurfData::addPoint.  Points in this data set "
	       << "have " << xsize << " dimensions and " << fsize
	       << " response values; point to be added has "
	       << sp.xSize() << " dimensions and " << sp.fSize()
	       << " response values. (Or gradient and Hessian sizes don't " 
	       << "match.)" << endl;
      throw bad_surf_data(errormsg.str());
    }
  }
  SurfPointSet::iterator iter;
  // This should be a safe const cast.  All that's happening is a check
  // to see if another data point at the same location in the space has already
  // been added.
  iter = orderedPoints.find(const_cast<SurfPoint*>(&sp));
  if (iter == orderedPoints.end()) {
    // This SurfPoint is not already in the data set.  Add it.
    points.push_back(new SurfPoint(sp));
    orderedPoints.insert(points[points.size()-1]);
    mapping.push_back(points.size()-1);
  } else {
    // Another SurfPoint in this SurfData object has the same location and
    // may have different response value(s).  Replace the old point with 
    // this new one.
    SurfPoint* spPtr = *iter;
    *spPtr = sp;
  }
}

/// Add a new response variable to each point. Return the index of the new 
/// variable.
unsigned SurfData::addResponse(const vector<double>& newValues, 
  string label)
{
  unsigned new_index;
  ostringstream errormsg;
  if (points.empty()) {
    throw bad_surf_data(
             "Cannot add response because there are no data points"
          );
  } else if (points.size() != mapping.size()) {
    errormsg << "Cannot add response because physical set size is different "
	     << "than logical set size.\nBefore adding another response, "
             << "clear \"excluded points\" or create a new data set by using " 
	     << "the SurfData::copyActive method." << endl;
    throw bad_surf_data(errormsg.str());
  } else if (newValues.size() != points.size()) {
    errormsg << "Cannot add another response: the number of new response"
             << " values does not match the size of the physical data set." 
             << endl;
    throw bad_surf_data(errormsg.str());
  } else {
    new_index = points[mapping[0]]->addResponse(newValues[0]);
    fsize++;
    for (unsigned i = 1; i < points.size(); i++) {
      new_index = points[mapping[i]]->addResponse(newValues[i]);
      assert(new_index == fsize - 1);
    }
  }
  if (label != "") {
    fLabels.push_back(label);
  } else {
    ostringstream labelos;
    labelos << "f" << new_index ;
    fLabels.push_back(labelos.str());
  }
  return new_index;
}

void SurfData::setConstraintPoint(const SurfPoint& sp)
{
  // handle the case of this being the first point in the data set
  if (points.empty()) { 
    xsize = sp.xSize();
    fsize = sp.fSize();
    gradsize = sp.fGradientsSize();
    hesssize = sp.fHessiansSize();
    if (xLabels.empty()) { 
      defaultLabels(); 
    } 
  } else { 
    if (sp.xSize() != xsize || sp.fSize() != fsize || 
	sp.fGradientsSize() != gradsize || sp.fHessiansSize() != hesssize) {
      ostringstream errormsg;
      errormsg << "Error in SurfData::setConstraintPoint.  Points in this data set "
               << "have " << xsize << " dimensions and " << fsize
               << " response values; point to be added has "
               << sp.xSize() << " dimensions and " << sp.fSize() 
               << " response values. (Or gradient and Hessian sizes don't " 
	       << "match.)" << endl;
      throw bad_surf_data(errormsg.str()); 
    } 
  } 
  constraintPoint = sp;
}

/// Specify which points should be skipped.  This can be used when only a 
/// subset of the SurfPoints should be used for some computation.
void SurfData::setExcludedPoints(const set<unsigned>& excluded_points)
{
  if (excluded_points.size() > points.size()) {
    throw bad_surf_data(
      "Size of set of excluded points exceeds size of SurfPoint set"
    );
  } else if (excluded_points.empty()) {
    defaultMapping();
    this->excludedPoints.clear();
  } else {
    // The size of the logical data set is the size of the physical
    // data set less the number of excluded points    
    mapping.resize(points.size() - excluded_points.size());
    unsigned mappingIndex = 0;
    unsigned sdIndex = 0;
    // map the valid indices to the physical points in points
    while (sdIndex < points.size()) {
      if (excluded_points.find(sdIndex) == excluded_points.end()) {
        mapping[mappingIndex++] = sdIndex;
      }
      sdIndex++;
    }
    this->excludedPoints = excluded_points;
    assert(mappingIndex == mapping.size());
  }
}

/// For use with copy constructor and assignment operator-- creates a list of
/// pointers to the points in the data set which is used to check for 
/// duplication when other points are added in the future
void SurfData::buildOrderedPoints()
{
  orderedPoints.clear();
  for (unsigned i = 0; i < points.size(); i++) {
    orderedPoints.insert(points[i]);
  }
}

/// Maps all indices to themselves in the mapping data member
void SurfData::defaultMapping()
{
  mapping.resize(points.size());
  for (unsigned i = 0; i < points.size(); i++) {
    mapping[i] = i;
  }
}
   
/// Set the labels for the predictor variables
void SurfData::setXLabels(const vector<string>& labels)
{
  if (labels.size() != xsize) {
    throw string("Dim mismatch in SurfData::setXLabels");
  }
  xLabels = labels;
}

/// Set the labels for the response variables
void SurfData::setFLabels(const vector<string>& labels)
{
  if (labels.size() != fsize) {
    throw string("Dim mismatch in SurfData::setFLabels");
  }
  fLabels = labels;
}

/// Set the label for a single response variable
void SurfData::setFLabel(unsigned index, const string& response_name)
{
  if (index >= fsize) {
    throw string("Dim mismatch in SurfData::setFLabel");
  }
  fLabels[index] = response_name;
}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write a set of SurfPoints to a file.  Opens the file and calls
/// other version of write.  All files include dimensions, but some
/// contain additional header and label information.
void SurfData::write(const string& filename) const
{
  if (mapping.empty()) {
    ostringstream errormsg;
    errormsg << "Cannot write SurfData object to stream."
	     << "  No active data points." << endl;
    throw bad_surf_data(errormsg.str());
  }
  bool binary = hasBinaryFileExtension(filename);
  ofstream outfile(filename.c_str(), 
    (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    writeBinary(outfile);
  } else {
    bool write_header = false;
    // Write the header and label info for .spd, not for .dat
    bool metadata = surfpack::hasExtension(filename,".spd");
    writeText(outfile, write_header, metadata);
  }
  outfile.close();
}

/// Read a set of SurfPoints from a file.  Opens file and calls other version.
void SurfData::read(const string& filename)
{
  // Open file in binary or text mode based on filename extension (.bspd or .spd)
  bool binary = hasBinaryFileExtension(filename);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    readBinary(infile);
  } else {
    readText(infile);
  }
  // Object may have already been created
  infile.close();
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeBinary(ostream& os) const 
{
  unsigned s = mapping.size();
  os.write((char*)&s,sizeof(s));
  os.write((char*)&xsize,sizeof(xsize));
  os.write((char*)&fsize,sizeof(fsize));
  os.write((char*)&gradsize,sizeof(gradsize));
  os.write((char*)&hesssize,sizeof(hesssize));
  ///\todo accumulate all the writes into one big chunk of data
  /// and write it out all at once.
  for (unsigned i = 0; i < mapping.size(); i++) {
    points[mapping[i]]->writeBinary(os);
  }
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os, 
			 bool write_header, bool write_labels) const
{
    if (write_header) {
      os << mapping.size() << endl
         << xsize << endl 
         << fsize << endl
         << gradsize << endl
         << hesssize << endl;
    }
    if (write_labels) {
      os << '%';
      for (unsigned i = 0; i < xLabels.size(); i++) {
        int correction = i ? 0 : 1; // the '%' takes up one space
        os << setw(surfpack::field_width - correction) << xLabels[i];
      }
      for (unsigned i = 0; i < fLabels.size(); i++) {
        os << setw(surfpack::field_width) << fLabels[i];
      }
      os << endl;
    }
    for (unsigned i = 0; i < mapping.size(); i++) {
      //if (!write_header) {
      //  os << setw((int)log10((double)(mapping.size())+2)) << (i+1);
      //}
      points[mapping[i]]->writeText(os);
    }
}

/// Read a set of SurfPoints from an input stream
void SurfData::readBinary(istream& is) 
{
  ///\todo Eliminate need for SurfPoint constructor that takes istream
  ///Instead, allocate a big chunk of memory, read all of the data in
  /// at the same time, then iterate through the block of data and
  /// create SurfPoints using the constructor that takes x and f.
  unsigned n_points_read = 0;
  unsigned size;
  try {
    cleanup();
    is.read((char*)&size,sizeof(size));
    is.read((char*)&xsize,sizeof(xsize));
    is.read((char*)&fsize,sizeof(fsize));
    is.read((char*)&gradsize,sizeof(gradsize));
    is.read((char*)&hesssize,sizeof(hesssize));
    points.clear();
    for (n_points_read = 0; n_points_read < size; n_points_read++) {
      // Throw an exception if we hit the end-of-file before we've
      // read the number of points that were supposed to be there.
      surfpack::checkForEOF(is);
      this->addPoint(SurfPoint(is, xsize, fsize, gradsize, hesssize));  
    }
    defaultMapping();
  } catch(surfpack::io_exception&) {
    cerr << "Expected: " << size << " points.  "
         << "Read: " << n_points_read << " points." << endl;
    throw;
  } 
}

unsigned SurfData::readHeaderInfo(istream& is)
{
  string single_line;

  getline(is,single_line);
  istringstream streamline(single_line);
  unsigned declared_size;
  streamline >> declared_size;

  getline(is,single_line);
  streamline.str(single_line); streamline.clear();
  streamline >> xsize;

  getline(is,single_line);
  streamline.str(single_line); streamline.clear();
  streamline >> fsize;

  getline(is,single_line);
  streamline.str(single_line); streamline.clear();
  streamline >> gradsize;

  getline(is,single_line);
  streamline.str(single_line); streamline.clear();
  streamline >> hesssize;
  
  return declared_size;
}

/// Read a set of SurfPoints from a text input stream
void SurfData::readText(istream& is, bool read_header, unsigned skip_columns) 
{
  unsigned n_points_read = 0;
  unsigned declared_size = 0;
  string single_line;
  try {
    cleanup();
    points.clear();
    if (read_header) declared_size = readHeaderInfo(is);

    getline(is,single_line);
    istringstream streamline(single_line);
    if (!readLabelsIfPresent(single_line)) {
      if (single_line != "" && single_line != "\n" && single_line[0] != '%') {
        this->addPoint(SurfPoint(single_line, xsize, 
				 fsize, gradsize, hesssize, skip_columns));
        n_points_read = 1;
      }
    }
    while (!is.eof()) {
      // Throw an exception if we hit the end-of-file before we've
      // read the number of points that were supposed to be there.
      //surfpack::checkForEOF(is);
      getline(is,single_line);
      if (single_line[0] == '%' || single_line == "") continue;
      // False for last argument signals a text read
      this->addPoint(SurfPoint(single_line, xsize, 
			       fsize, gradsize, hesssize, skip_columns));  
      n_points_read++;
    }
    defaultMapping();
  } catch(surfpack::io_exception& exception) {
    cerr << exception.what() << endl;
    throw;
  } 
  if (read_header && n_points_read != declared_size) {
    ostringstream errormsg;
    errormsg << "Expected: " << declared_size << " points.  "
         << "Read: " << n_points_read << " points." << endl;
    throw surfpack::io_exception(errormsg.str());
  }
}

// Print set of data points to a stream. 
ostream& operator<<(ostream& os, const SurfData& sd) 
{ 
  sd.writeText(os); 
  return os;
}

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

/// Returns true if file has .bspd extension, false if it has .spd extension. 
/// Otherwise, an exception is thrown.
bool SurfData::hasBinaryFileExtension(const string& filename) const
{
  if (surfpack::hasExtension(filename,".bspd")) {
    return true;
  } else if (surfpack::hasExtension(filename,".spd")) {
    return false;
  } else if (surfpack::hasExtension(filename,".dat")) {
    return false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .bspd, or .spd"
    );
  }
}

/// Set x vars labels to x0 x1, etc.; resp. vars to f0 f1, etc.
void SurfData::defaultLabels()
{
  xLabels.resize(xsize);
  for (unsigned i = 0; i < xsize; i++) {
    ostringstream os;
    os << "x" << i ;
    xLabels[i] = os.str();
  }
  fLabels.resize(fsize);
  for (unsigned i = 0; i < fsize; i++) {
    ostringstream os;
    os << "f" << i ;
    fLabels[i] = os.str();
  }
}

bool SurfData::readLabelsIfPresent(string single_line)
{
  if (single_line[0] != '%') {
    defaultLabels();
    return false;
  } else {
    single_line[0] = ' ';
    string label;
    xLabels.resize(xsize);
    istringstream is(single_line);
    for (unsigned i = 0; i < xsize; i++) {
      is >> xLabels[i];
      if (xLabels[i] == "") { 
        // not enough heading names 
        // line of column headings.  Use the default headings and return.
        defaultLabels();
        return false;
      }
    } // predictor variable labels
    fLabels.resize(fsize);
    for (unsigned i = 0; i < fsize; i++) {
      is >> fLabels[i];
      if (fLabels[i] == "") { 
        // not enough heading names 
        // line of column headings.  Use the default headings and return.
        defaultLabels();
        return false;
      }
    } // response variable labels
  } // custom labels
  return true;
}
// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

// Throw an exception if there are any mismatches in the number of
// dimensions or number of response values among points in the data set  
void SurfData::sanityCheck() const
{
  if (!points.empty()) {
    unsigned dimensionality = points[0]->xSize();
    unsigned numResponses = points[0]->fSize();
    unsigned numGrad = points[0]->fGradientsSize();
    unsigned numHess = points[0]->fHessiansSize();
    for (unsigned i = 1; i < points.size(); i++) {
      if (points[i]->xSize() != dimensionality ||
          points[i]->fSize() != numResponses || 
          points[i]->fGradientsSize() != numGrad || 
          points[i]->fHessiansSize() != numHess ) {
        ostringstream errormsg;
        errormsg << "Error in SurfData::sanityCheck." << endl
		 << "Point 0 has " << dimensionality << " dimensions "
                 << "and " << numResponses << " response values, " << endl
                 << "but point " << i << " has " << points[i]->xSize()
 		 << " dimensions and " << points[i]->fSize() << "response "
 		 << " values. (Or gradient and Hessian sizes are wrong.)";
	throw bad_surf_data(errormsg.str());
      } // if mismatch
    } // for each point
  } // if !empty
}

/// Check that the index falls within acceptable boundaries (i.e., is
/// less than mapping.size()
void SurfData::checkRangeNumPoints(const string& header, unsigned index) const
{
  if (index >= mapping.size()) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (mapping.empty()) {
      errormsg << "Index " << index << " specified, but there are zero points "
	       << "in the logical data set."
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << mapping.size() - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}

/// Make sure an index falls within acceptable boundaries (i.e., index is less
/// than fsize)
void SurfData::checkRangeNumResponses(const string& header, 
  unsigned index) const
{
  if (index >= fsize) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (fsize == 0) {
      errormsg << "Index " << index 
	       << " specified, but there are zero response"
	       << "values."
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << fsize - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}

VecVecDbl SurfData::asVecVecDbl(const SurfData& data)
{
  VecVecDbl result(data.size());
  for (unsigned i = 0; i < data.size(); i++) {
    result[i].resize(data.xSize());
    for (unsigned j = 0; j < data.xSize(); j++) {
      result[i][j] = data(i,j);
    }
  }
  return result;
}
