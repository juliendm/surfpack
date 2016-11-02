/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_SYSTEM_HEADERS_H
#define SURFPACK_SYSTEM_HEADERS_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  // HAVE_CONFIG_H is STILL set in Dakota/src (EVEN IN THE CMAKE BUILD!) so
  // use a "disable config header" conditional to help manage the transition
  #include "surfpack_config.h"
#endif // HAVE_CONFIG_H

#include <algorithm>
#include <limits>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>
#include <exception>

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

// crucial that these be included before any calls to BOOST_CLASS_EXPORT
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp> 
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/access.hpp>

// serialization for STL types
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

#endif

typedef std::vector<double> VecDbl;
typedef std::vector<double>::const_iterator VecDblIt;
typedef std::vector<unsigned> VecUns;
typedef std::vector<unsigned>::const_iterator VecUnsIt;
typedef std::vector< std::vector< unsigned > > VecVecUns;
typedef std::vector< std::vector< double > > VecVecDbl;
typedef std::vector< std::string > VecStr;
typedef std::set< unsigned > SetUns;

typedef std::pair< std::string, std::string > ModelParam;
typedef std::map< std::string, std::string> ParamMap;
typedef std::pair< std::string, ParamMap> Command;

#endif // SURFPACK_SYSTEM_HEADERS_H

