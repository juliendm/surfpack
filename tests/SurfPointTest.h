/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef SURFPOINTTEST_H
#define SURFPOINTTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "SurfPoint.h"
#include "SurfData.h"

class SurfPointTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfPointTest );
  CPPUNIT_TEST( testSurfPointPtrLessThan );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testConstructorXSpecified );
  CPPUNIT_TEST( testConstructorXSpecifiedPlusOneF );
  CPPUNIT_TEST( testConstructorXSpecifiedFVector );
  CPPUNIT_TEST( testConstructorFromIStreamBinary );
  CPPUNIT_TEST( testConstructorFromIStreamText );
  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST_EXCEPTION( testConstructorBadXSize, SurfPoint::null_point );
  CPPUNIT_TEST( testOperatorAssignment );
  CPPUNIT_TEST( testOperatorAssignmentToSelf );
  CPPUNIT_TEST( testOperatorEquality );
  CPPUNIT_TEST( testOperatorInequality );
  CPPUNIT_TEST( testXSize );
  CPPUNIT_TEST( testFSize );
  CPPUNIT_TEST( testX );
  CPPUNIT_TEST( testFQuery );
  CPPUNIT_TEST_EXCEPTION( testFQueryBadIndex, std::range_error );
  CPPUNIT_TEST( testAddResponse );
  CPPUNIT_TEST( testFAssign );
  CPPUNIT_TEST_EXCEPTION( testFAssignBadIndex , std::range_error );
  CPPUNIT_TEST( testWriteBinary );
  CPPUNIT_TEST( testWriteText );
  CPPUNIT_TEST( testReadBinary );
  CPPUNIT_TEST( testReadText );
  CPPUNIT_TEST( testStreamInsertion );
  CPPUNIT_TEST( testResize );
  CPPUNIT_TEST( testSetX );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();

// Other
  void testSurfPointPtrLessThan();

// Constructors
  void testConstructor();
  void testConstructorXSpecified();
  void testConstructorXSpecifiedPlusOneF();
  void testConstructorXSpecifiedFVector();
  void testConstructorFromIStreamBinary();
  void testConstructorFromIStreamText();
  void testCopyConstructor();
  void testConstructorBadXSize();

// Overloaded operators
  void testOperatorAssignment();
  void testOperatorAssignmentToSelf();
  void testOperatorEquality();
  void testOperatorInequality();

// Queries
  void testXSize();
  void testFSize();
  void testX();
  void testFQuery();
  void testFQueryBadIndex();

// Commands
  void testAddResponse();
  void testFAssign();
  void testFAssignBadIndex();
  void testResize();
  void testSetX();

// I/O
  void testWriteBinary();
  void testWriteText();
  void testReadBinary();
  void testReadText();
  void testStreamInsertion();

private:
  std::vector<double> x1;
  std::vector<double> x2;
  std::vector<double> x3;
  std::vector<double> f1;
  std::vector<double> f2;

  SurfPoint* spPtr;
  SurfPoint* spPtr2;
 
  
};


#endif
