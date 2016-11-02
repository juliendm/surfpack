/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef SURFPACKCOMMONTEST_H
#define SURFPACKCOMMONTEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "SurfPoint.h"
#include "SurfData.h"

class SurfpackCommonTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfpackCommonTest );
CPPUNIT_TEST( testConstructor );
CPPUNIT_TEST( matrixMultiplyTest );
CPPUNIT_TEST( matrixVectorMultTest );
CPPUNIT_TEST( matrixVectorMultTransTest );
CPPUNIT_TEST( matrixReshapeTest );
CPPUNIT_TEST( matrixResizeCTest );
CPPUNIT_TEST( stringToVecUnsTest );
CPPUNIT_TEST( weightedAvgTest );
CPPUNIT_TEST( toString );
CPPUNIT_TEST( fromVec );
CPPUNIT_TEST( blockTests );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void testConstructor();
void matrixMultiplyTest();
void matrixVectorMultTest();
void matrixVectorMultTransTest();
void matrixReshapeTest();
void matrixResizeCTest();
void stringToVecUnsTest();
void weightedAvgTest();
void toString();
void fromVec();
void blockTests();
};

#endif
