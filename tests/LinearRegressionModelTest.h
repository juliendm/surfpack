/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef LINEAR_REGRESSION_MODEL_TEST_H 
#define LINEAR_REGRESSION_MODEL_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "SurfPoint.h"
#include "SurfData.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "LinearRegressionModel.h"

class LinearRegressionModelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LinearRegressionModelTest );
//CPPUNIT_TEST( constructorTest );
//CPPUNIT_TEST( unityBasisTest );
//CPPUNIT_TEST( singleLinearTest );
//CPPUNIT_TEST( singleQuadraticTest );
//CPPUNIT_TEST( lineEvalTest );
//CPPUNIT_TEST( quadratic2DTest );
//CPPUNIT_TEST( plotTest1 );
CPPUNIT_TEST( termPrinterTest );
CPPUNIT_TEST( createModelTest );
//CPPUNIT_TEST( FTest );
  CPPUNIT_TEST_SUITE_END();
public:
  AxesBounds* ab;
  SurfData* sd;
  SurfData* randsd;
  LRMBasisSet bs;
  void setUp();
  void tearDown();
void constructorTest();
void lineEvalTest();
void unityBasisTest();
void singleLinearTest();
void singleQuadraticTest();
void quadratic2DTest();
void plotTest1();
void termPrinterTest();
void createModelTest();
//void FTest();
};

#endif
