/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef UNITTESTS_H
#define UNITTESTS_H

#include <string>

const unsigned unsignedZero = 0;
const double doubleZero = 0.0;
const int intZero = 0;

void writePoint1Files();
void writePoint2Files();
void writeRastriginAndClaimsTooManyFiles();
void writeManyPtsFiles();
void writeOneDimQuadratic();
void writeUnknownSurfaceFile();
void writeSurfDataWithHeader();
void setOstreamFlags(std::ostream& os);
void initialize();
bool matches(double observed, double target, double margin = 1.0e-2);

#endif
