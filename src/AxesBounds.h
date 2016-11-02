/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef AXES_BOUNDS_H
#define AXES_BOUNDS_H

#include "surfpack_system_headers.h"

class SurfData;
class SurfPoint;

/// Concrete class used in conjunction with Surfpack command CreateSample. 
/// Minimum and maximum values are specified along each dimension in the 
/// parameter space.  A SurfData object can be created by requesting
/// simple Monte Carlo samples from hypercube defined by the boundaries.
/// Alternatively, a client may specify a number of points along each dimension
/// (in addition to the maximum and minimum) and then create a SurfData object 
/// from that hypergrid of points.
class AxesBounds
{
public:
  /// Values for one dimension For random samples, only bounds are used.  For
  /// hyper-gridding, the pts and interval fields are also used.  
  struct Axis {
      Axis() : min(0.0), max(0.0), minIsMax(true) {}
      Axis(double min_in, double max_in) : min(min_in), max(max_in), minIsMax(min==max) {}

      /// Minimum value along this dimension.
      double min;

      /// Maximum value along this dimension.
      double max;

      /// No variation along this dimension
      bool minIsMax;

      std::string asString() const;
  };

  /// Object is created from an existing list of Axis objects
  AxesBounds( std::vector<Axis> axes_in);

  /// Object is created from a string containing min/max pairs
  // along one or more axes, delimited by | 
  AxesBounds( std::string bounds);

  /// No special behavior
  ~AxesBounds();

  /// Return a hypergrid data set as a SurfData object.  The client is 
  /// responsible to deallocate the memory.
  SurfData* sampleGrid(const std::vector<unsigned>& grid_points, 
    const std::vector<std::string>& test_functions) const;
  SurfData* sampleGrid(const std::vector<unsigned>& grid_points) const; 

  /// Return a data set with size SurfPoints.  Parameter test_functions
  /// must contain the names of zero or more functions at which all the data
  /// points should be evaluated.  The test functions should reside in the
  /// surfpack namespace.  The client must deallocate the memory.
  SurfData* sampleMonteCarlo(unsigned size, 
    const std::vector<std::string>& test_functions) const;
  SurfData* sampleMonteCarlo(unsigned size) const; 

  /// Advance the counter used to iterate through dimensions for grid data.
  /// For example, if the client has requested a 10 x 10 grid, and the point
  /// data member currently holds (3,9), it will be advanced (not unlike an
  /// odometer) to (4,0).  This will signify that the next point to be added
  /// should be (axes[0].min+4*axes[0].interval, axes[1].min+0*axes[1].interval.
  void nextPoint(std::vector<unsigned>& point_odometer,
    const std::vector<unsigned>& grid_points) const;

  std::vector<double> computeIntervals(const std::vector<Axis>& axes, 
    const std::vector<unsigned>& grid_points) const;

  /// Return the information stores in m_axes as a text string
  std::string asString() const;

  /// Return the number of dimensions 
  unsigned size() const;

  /// Return a reference to the bounds for the dimension specified by index
  const Axis& operator[](unsigned index) const;

  /// Return a reference to the vector of Axis
  const std::vector<Axis>& axes() const;

  /// Return the smallest hypercube that contains all the data
  static AxesBounds boundingBox(const SurfData& sd);

protected:
  /// Parse out the boundary information from an input stream created from a 
  /// string or file
  void parseBounds(std::istream& is);

protected:
  /// The set of <minimum, maximum> specifications for each dimension
  std::vector<Axis> m_axes;
};
#endif
