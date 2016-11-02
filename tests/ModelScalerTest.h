/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef MODEL_SCALER_TEST_H 
#define MODEL_SCALER_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "SurfPoint.h"
#include "SurfData.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "ModelScaler.h"
class ModelScalerTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ModelScalerTest );
CPPUNIT_TEST( nonScaleTest );
CPPUNIT_TEST( NormalizingScalerDataTest );
CPPUNIT_TEST( NormalizingScalerModelTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();

void nonScaleTest();
void NormalizingScalerDataTest();
void NormalizingScalerModelTest();
};

#endif
