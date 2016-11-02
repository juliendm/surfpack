/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif


#ifndef RADIAL_BASIS_FUNCTION_TEST_H 
#define RADIAL_BASIS_FUNCTION_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "Conmin.h"
#include "SurfpackModel.h"

class RadialBasisFunctionTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( RadialBasisFunctionTest );
//  CPPUNIT_TEST( generalModelTest );
//CPPUNIT_TEST( hillValleyTest );
//CPPUNIT_TEST( eTest );
//CPPUNIT_TEST( partitionTest );
//CPPUNIT_TEST( centroidTest );
//CPPUNIT_TEST( centroidTest2 );
//CPPUNIT_TEST( updateCentroidTest );
//CPPUNIT_TEST( cvtTest );
CPPUNIT_TEST( createTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
void generalModelTest();
void hillValleyTest();
void eTest();
void partitionTest();
void centroidTest();
void centroidTest2();
void updateCentroidTest();
void cvtTest();
void createTest();
};

#endif
