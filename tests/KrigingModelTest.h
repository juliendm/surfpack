/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef KRIGING_MODEL_TEST_H 
#define KRIGING_MODEL_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

class KrigingModelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( KrigingModelTest );
CPPUNIT_TEST( simpleTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
void simpleTest();
};

#endif
