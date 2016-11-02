/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::istream;
using std::istringstream;
using std::multiplies;
using std::string;
using std::vector;
using surfpack::dbg;
const int axdbg = 0;

  
std::string AxesBounds::Axis::asString() const
{ 
  std::ostringstream os;
  os << min;
  if (!minIsMax) os << " " << max;
  return os.str();
}
  
unsigned AxesBounds::size() const
{
  return m_axes.size();
}

const AxesBounds::Axis& AxesBounds::operator[](unsigned index) const
{
  assert(index < size());
  return m_axes[index];
}

const std::vector<AxesBounds::Axis>& AxesBounds::axes() const
{
  return m_axes;
}

/// Object is created from an existing list of Axis objects
AxesBounds::AxesBounds(vector<AxesBounds::Axis> axes_in) 
  : m_axes(axes_in)
{

}

/// Object is created from a text file
AxesBounds::AxesBounds(string bounds) 
{
    // Consider a string with the following data:
    // 0 1 | 2 | -1 1

    // Each pair or singleton corresponds to one variable
    // If there are two values, they represent the minimum
    // and maximum value for that variable.  If there is only one value
    // that means that the variable is fixed (fixed variable = oxymoron?).
    // White space is actually ignored.  The vertical bars between the values
    // for different dimensions are required.

    istringstream is(bounds + " ");
    parseBounds(is);
}    

/// Read (min,max) pairs or a fixed value for each dimension.
/// Data for different dimensions should be separated by '|' 
void AxesBounds::parseBounds(istream& is) 
{
  m_axes.push_back(Axis());
  string token;
  while (!is.eof()) {
    // read in min for current dimension
    is >> m_axes.back().min;
    // read in next token -- it should be either a '|' or the max for this dim
    is >> token;
    dbg(axdbg) << "Token read; " << token << " eof?: " << is.eof() << "\n";
    if (is.eof()) break;
    if (token == "|") {
      m_axes.back().max = m_axes.back().min;
      m_axes.push_back(Axis());
      continue; // There is no 'max' for this dim; skip to next one
    } else {
      m_axes.back().max = std::atof(token.c_str());
      m_axes.back().minIsMax = false;
      // now the next token should be a '|' or eof
      is >> token;
    dbg(axdbg) << "Token read; " << token << " eof?: " << is.eof() << "\n";
      if (is.eof()) break;
      if (token != "|") {
        cerr << "Expected |" << endl;
 	std::exit(1);
      }
      m_axes.push_back(Axis());
    }
  }

  if (axdbg) { // debug output
    cout << "Axes values parsed" << endl;
    for (unsigned i = 0; i < m_axes.size(); i++) {
      cout << m_axes[i].min;
      if (!m_axes[i].minIsMax) {
        cout << " " << m_axes[i].max;
      }
      cout << endl;
    } 
    cout << "dims: " << m_axes.size() << endl;
  } // debug output
      
}

/// Based on the ranges for each variable contained in m_axes and the number
/// of grid points requested per dimension (grid_points), computes and returns t/// the interval between two values on each dimension.
vector<double> AxesBounds::computeIntervals(const vector<Axis>& axes, 
  const vector<unsigned>& grid_points) const
{
  assert(axes.size() == grid_points.size());
  vector<double> intervals(axes.size());
  for (unsigned i = 0; i < grid_points.size(); i++) {
    // Treating one as special case avoids div. by zero below
    if (grid_points[i] == 1) {
      intervals[i] = 0.0;
    } else {
      dbg(axdbg) << "i " << i << " min/max: " << axes[i].min << " " << axes[i].max << " gp: " << grid_points[i] << " int: ";
      intervals[i] = (axes[i].max - axes[i].min)/(grid_points[i] - 1);
      dbg(axdbg) << intervals[i] << "\n";
    }
  }
  return intervals;
}     

/// No special behavior
AxesBounds::~AxesBounds()
{

}


/// Advance the counter used to iterate through dimensions for grid
/// data.  For example, if the client has requested a 10 x 4 grid, and
/// the point_odometer currently holds (7,3), it will be advanced
/// (like an odometer) to (8,0).  This will signify that the next
/// point to be added should be (m_axes[0].min+8*intervals[0],
/// m_axes[1].min+0*intervals[1]).
///
/// NOTE: On the last call the odometer will increase beyond the 0th
/// grid_points bound; the caller is responsible for not requesting
/// too many points and exceeding the bounds.
void AxesBounds::nextPoint(vector<unsigned>& point_odometer, 
  const vector<unsigned>& grid_points) const
{
  // Scan across the "odometer" reading to find the first value that
  // does not need to roll over.  For example, if the user requested a
  // 5 x 5 x 5 x 5 grid and point_odometer holds (3,2,4,4), then the
  // next value should be (3,3,0,0), so the rightmost two values need
  // to roll over, and cur_dim should point to the 2's position.  Skip
  // over dimensions with one grid point (fixed).
  int cur_dim = m_axes.size()-1;
  while (cur_dim > 0 && 
	 (grid_points[cur_dim] == 1 || 
	  point_odometer[cur_dim] == (grid_points[cur_dim] - 1))
	 ) {
    cur_dim--;
  }

  // If the odometer isn't maxed out, increase the digit at cur_dim by one,
  // and then zero out everything to the right.
  point_odometer[cur_dim]++;
  cur_dim++;
  while(cur_dim < m_axes.size()) {
    point_odometer[cur_dim] = 0;
    cur_dim++;
  }
}

/// Return a hypergrid data set as a SurfData object.  The client is 
/// responsible to deallocate the memory.
SurfData* AxesBounds::sampleGrid(const vector<unsigned>& grid_points) const
{
  return sampleGrid(grid_points,vector<string>());
}

/// Return a hypergrid data set as a SurfData object, evaluating
/// requested test_function, if any.
SurfData* AxesBounds::sampleGrid(const vector<unsigned>& grid_points, 
				 const vector<string>& test_functions) const
{
  vector<unsigned> point_odometer(grid_points.size(),0);
  vector<double> surfptx(m_axes.size());
  vector<SurfPoint> sps;
  vector<double> intervals = computeIntervals(m_axes,grid_points);
  unsigned npts = accumulate(grid_points.begin(),grid_points.end(),
	1, multiplies<unsigned>());
  for (int i = 0; i < npts; i++) {
      for (int j = 0; j < m_axes.size(); j++) {
          surfptx[j] = m_axes[j].min + intervals[j]*point_odometer[j];
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < test_functions.size(); k++) {
        sp.addResponse(surfpack::testFunction(test_functions[k], sp.X()));
      }
      sps.push_back(sp);
      nextPoint(point_odometer, grid_points);
  }
  SurfData* sd = new SurfData(sps);
  if (!test_functions.empty()) sd->setFLabels(test_functions);
  return sd;
}

/// Return a data set with size SurfPoints.  Parameter test_functions
/// must contain the names of zero or more functions at which all the data
/// points should be evaluated.  The test functions should reside in the
/// surfpack namespace.
SurfData* AxesBounds::sampleMonteCarlo(unsigned size) const
{
  return sampleMonteCarlo(size,vector<string>());
}
SurfData* AxesBounds::sampleMonteCarlo(unsigned size, 
  const vector<string>& test_functions) const
{
  vector<double> surfptx(m_axes.size());
  vector<SurfPoint> sps;
  for (unsigned i = 0; i < size; i++) {
      for (unsigned j = 0; j < m_axes.size(); j++) {
	surfptx[j] = (m_axes[j].max - m_axes[j].min) * 
          (surfpack::shared_rng().rand())+ m_axes[j].min;
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < test_functions.size(); k++) {
        sp.addResponse(surfpack::testFunction(test_functions[k], sp.X()));
      }
      sps.push_back(sp);
  }
  return new SurfData(sps);
}

std::string AxesBounds::asString() const
{
  std::ostringstream os;
  for (unsigned i = 0; i < m_axes.size(); i++) {
    os << m_axes[i].asString() << "\n";
  }
  return os.str();
}

AxesBounds AxesBounds::boundingBox(const SurfData& sd)
{
  assert(sd.size());
  assert(sd.xSize());
  vector<Axis> axes(sd.xSize()); 
  for (unsigned i = 0; i < axes.size(); i++) {
    axes[i].min = std::numeric_limits<double>::max();
    axes[i].max = -std::numeric_limits<double>::max();
  }
  for (unsigned pt = 0; pt < sd.size(); pt++) {
    for (unsigned dim = 0; dim < sd.xSize(); dim++) {
      if (sd(pt,dim) < axes[dim].min) {
        axes[dim].min = sd(pt,dim);
      }
      if (sd(pt,dim) > axes[dim].max) {
        axes[dim].max = sd(pt,dim);
      }
    }
  }
  for (unsigned i = 0; i < axes.size(); i++) {
    if (axes[i].min != axes[i].max) axes[i].minIsMax = false;
  }
  return AxesBounds(axes);
}
