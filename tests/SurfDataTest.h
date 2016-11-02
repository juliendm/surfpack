/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef SURFDATATEST_H
#define SURFDATATEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "surfpack.h"
#include "SurfData.h"

class SurfDataTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfDataTest );
  CPPUNIT_TEST( testConstructorVectorPoints );
  CPPUNIT_TEST( testConstructorVectorPointsEmpty);
  CPPUNIT_TEST_EXCEPTION( testConstructorVectorPointsMismatchXs, 
    SurfData::bad_surf_data); 
  CPPUNIT_TEST_EXCEPTION( testConstructorVectorPointsMismatchFs, 
    SurfData::bad_surf_data);
  CPPUNIT_TEST( testConstructorFilenameText );
  CPPUNIT_TEST( testConstructorFilenameBinary );
  CPPUNIT_TEST( testConstructorIStreamText );
  CPPUNIT_TEST( testConstructorIStreamBinary );
  CPPUNIT_TEST( testCopyConstructorSimple );
  CPPUNIT_TEST( testCopyConstructorComplex );
  CPPUNIT_TEST( testCopyActive );
  CPPUNIT_TEST( testCopyActiveEmpty );
  CPPUNIT_TEST( testAssignment );
  CPPUNIT_TEST( testAssignmentToSelf );
  CPPUNIT_TEST( testOperatorEquality );
  CPPUNIT_TEST( testOperatorInequality );
  CPPUNIT_TEST( testOperatorIndexing );
  CPPUNIT_TEST( testOperatorIndexingScaled );
  CPPUNIT_TEST_EXCEPTION( testOperatorIndexingBadIndex, std::range_error );
  CPPUNIT_TEST_EXCEPTION(testOperatorIndexingAnotherBadIndex, std::range_error);
  CPPUNIT_TEST( testSize );
  CPPUNIT_TEST( testEmpty );
  CPPUNIT_TEST( testXSize );
  CPPUNIT_TEST( testFSize );
  CPPUNIT_TEST( testGetExcludedPoints );
  CPPUNIT_TEST( testGetResponse );
  CPPUNIT_TEST( testGetDefaultIndex );
  CPPUNIT_TEST( testHasBinaryExtension);
  CPPUNIT_TEST( testHasTextExtension);
  CPPUNIT_TEST( testSetDefaultIndex);
  CPPUNIT_TEST_EXCEPTION( testSetDefaultIndexBadIndex, std::range_error);
  CPPUNIT_TEST_EXCEPTION( testSetDefaultIndexNoResponses, std::range_error);
  CPPUNIT_TEST( testSetResponse );
  CPPUNIT_TEST_EXCEPTION( testSetResponseBadIndex, std::range_error );
  CPPUNIT_TEST_EXCEPTION( testSetResponseAnotherBadIndex, 
    std::range_error );
  CPPUNIT_TEST( testAddPoint );
  CPPUNIT_TEST_EXCEPTION( testAddPointBadDimension, 
    SurfData::bad_surf_data );
  CPPUNIT_TEST_EXCEPTION( testAddPointBadNumResponses, 
    SurfData::bad_surf_data );
  CPPUNIT_TEST( testAddPointToEmptySet );
  CPPUNIT_TEST( testAddResponse );
  CPPUNIT_TEST_EXCEPTION( testAddResponseToEmptySet,
    SurfData::bad_surf_data );
  CPPUNIT_TEST_EXCEPTION( testAddResponseWithSkipped,
    SurfData::bad_surf_data );
  CPPUNIT_TEST_EXCEPTION( testAddResponseWrongNumber,
    SurfData::bad_surf_data );
  CPPUNIT_TEST( testSetExcludedPoints);
  CPPUNIT_TEST( testSetExcludedPointsToNone);
  CPPUNIT_TEST_EXCEPTION( testSetExcludedPointsTooMany,
    SurfData::bad_surf_data );
  CPPUNIT_TEST( testDuplicatePoint );
  CPPUNIT_TEST( testSetScalerNull );
  CPPUNIT_TEST( testSetScalerNotNull );
  CPPUNIT_TEST( testIsScaled );
  CPPUNIT_TEST( testWriteBinary );
  CPPUNIT_TEST( testWriteText );
  CPPUNIT_TEST_EXCEPTION( testWriteNoPoints,
    SurfData::bad_surf_data );
  CPPUNIT_TEST_EXCEPTION( testWriteNoFile, surfpack::file_open_failure);
  CPPUNIT_TEST_EXCEPTION( testWriteBadFileExtension,
    surfpack::io_exception);
  CPPUNIT_TEST_EXCEPTION( testReadNoFile, surfpack::file_open_failure);
  CPPUNIT_TEST_EXCEPTION( testReadBinaryFileTooShort,
    surfpack::io_exception);
  CPPUNIT_TEST_EXCEPTION( testReadTextFileTooShort,
    surfpack::io_exception);
  CPPUNIT_TEST_EXCEPTION( testBadSanityCheck,
    SurfData::bad_surf_data);
  CPPUNIT_TEST( testStreamInsertion );
CPPUNIT_TEST( columnHeaderTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void testConstructorVectorPoints();
  void testConstructorVectorPointsEmpty();
  void testConstructorVectorPointsMismatchXs();
  void testConstructorVectorPointsMismatchFs();
  void testConstructorFilenameText();
  void testConstructorFilenameBinary();
  void testConstructorIStreamText();
  void testConstructorIStreamBinary();
  void testCopyConstructorSimple();
  void testCopyConstructorComplex();
// Destructor not explicitly tested
// init not explicitly tested
  void testCopyActive();
  void testCopyActiveEmpty();
// copyBlockData not explicitly tested
// cleanup not explicitly tested
  void testAssignment();
  void testAssignmentToSelf();
  void testOperatorEquality();
  void testOperatorInequality();
  void testOperatorIndexing();
  void testOperatorIndexingScaled();
  void testOperatorIndexingBadIndex();
  void testOperatorIndexingAnotherBadIndex();
  void testSize();
  void testEmpty();
  void testXSize();
  void testFSize();
  void testGetExcludedPoints();
  void testGetResponse();
  void testGetDefaultIndex();
  void testHasBinaryExtension();
  void testHasTextExtension();
  void testSetDefaultIndex();
  void testSetDefaultIndexBadIndex();
  void testSetDefaultIndexNoResponses();
  void testSetResponse();
  void testSetResponseBadIndex();
  void testSetResponseAnotherBadIndex();
  void testAddPoint();
  void testAddPointBadDimension();
  void testAddPointBadNumResponses();
  void testAddPointToEmptySet();
  void testAddResponse();
  void testAddResponseToEmptySet();
  void testAddResponseWithSkipped();
  void testAddResponseWrongNumber();
  void testSetExcludedPoints();
  void testSetExcludedPointsToNone();
  void testSetExcludedPointsTooMany();
  void testDuplicatePoint();
  void testSetScalerNull();
  void testSetScalerNotNull();
  void testIsScaled();
// defaultMapping not explicitly tested
// validateXMatrix not explicitly tested
// validateYVector not explicitly tested
  void testWriteBinary();
  void testWriteText();
  void testWriteNoPoints();
  void testWriteNoFile();
  void testWriteBadFileExtension();
  void testReadNoFile();
  void testReadTextFileTooShort();
  void testReadBinaryFileTooShort();
  void testBadSanityCheck();
  void testStreamInsertion();

private:
  SurfData* sdPtr1;
  SurfData* sdPtr2;
  std::vector<SurfPoint> surfpoints;
  std::set<unsigned> skipAllPoints;
  std::set<unsigned> skipPoints;
void columnHeaderTest();
};

#endif
