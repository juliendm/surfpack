/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef MOVING_LEAST_SQUARES_TEST_H 
#define MOVING_LEAST_SQUARES_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "Conmin.h"
#include "SurfpackModel.h"

class MovingLeastSquaresTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MovingLeastSquaresTest );
//CPPUNIT_TEST( generalModelTest );
CPPUNIT_TEST( sineCurve );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
void generalModelTest();
void sineCurve();
};

#endif
