/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef SURFPACK_MODEL_TEST_H 
#define SURFPACK_MODEL_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

#include "SurfPoint.h"
#include "SurfData.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackModel.h"
class SurfpackModelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfpackModelTest );
CPPUNIT_TEST( plotTest1 );
//CPPUNIT_TEST( graphicalDerivTest );
//CPPUNIT_TEST( modelSampleTest );
CPPUNIT_TEST( manualANNTest );
  CPPUNIT_TEST_SUITE_END();
public:
  AxesBounds* ab;
  SurfData* sd;
  SurfData* randsd;
  LRMBasisSet bs;
  void setUp();
  void tearDown();

static void generateThePlot(const std::string& plotname, const std::string& datafilename, int x, int y1, int y2);
static void generalDerivativeTest(const SurfpackModel& model, const AxesBounds& ab);
void plotTest1();
void graphicalDerivTest();
void modelSample(const SurfpackModel& model);
void modelSampleTest();
void manualANNTest();
};

#endif
