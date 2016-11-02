/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>

#include "LinearRegressionModel.h"
#include "ModelScalerTest.h"
#include "ModelScaler.h"
#include "KrigingModel.h"
#include "DirectANNModel.h"
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"
#include "SurfpackInterface.h"
#include "unittests.h"
#include "ModelFitness.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::vector;
using std::ostream_iterator;
using std::ostringstream;

CPPUNIT_TEST_SUITE_REGISTRATION( ModelScalerTest );


void ModelScalerTest::setUp()
{

}

void ModelScalerTest::tearDown()
{
}

void ModelScalerTest::nonScaleTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","11 11","sphere");  
  //AxesBounds* ab = new AxesBounds("-2 2 | -2 2");
  //SurfData* sd = SurfpackInterface::CreateSample(ab, VecUns(11,11));
  // TODO: sphere
  CPPUNIT_ASSERT(sd->size() == 121);
  VecDbl responses = sd->getResponses();
  delete sd;
  double max_response = *max_element(responses.begin(),responses.end());
  CPPUNIT_ASSERT(matches(max_response,8.0));
  double min_response = *min_element(responses.begin(),responses.end());
  CPPUNIT_ASSERT(matches(min_response,0.0));
}

void ModelScalerTest::NormalizingScalerDataTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","11 11","sphere");
  //AxesBounds* ab = new AxesBounds("-2 2 | -2 2");
  //SurfData* sd = SurfpackInterface::CreateSample(ab, VecUns(11,11));
  CPPUNIT_ASSERT(sd->size() == 121);
  VecDbl pt(2,0.0);
  ModelScaler* ms = NormalizingScaler::Create(*sd);
  // Min is -2 and range is 4
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ms)->scalers[0].offset,-2.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ms)->scalers[1].offset,-2.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ms)->scalers[0].scaleFactor,4.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ms)->scalers[1].scaleFactor,4.0));
  ScaledSurfData ssd(*ms,*sd);
  VecDbl responses = ssd.getResponses();
  delete ms;
  delete sd;
  double max_response = *max_element(responses.begin(),responses.end());
  CPPUNIT_ASSERT(matches(max_response,1.0));
  double min_response = *min_element(responses.begin(),responses.end());
  CPPUNIT_ASSERT(matches(min_response,0.0));
}

void ModelScalerTest::NormalizingScalerModelTest()
{
  AxesBounds* ab = new AxesBounds("-10 10 | -10 10");
  SurfData* sd = SurfpackInterface::CreateSample(ab, VecUns(11,11));
  //  SurfData* sd = SurfpackInterface::CreateSample("-10 10 | -10 10","11 11","sphere");
  CPPUNIT_ASSERT(sd->size() == 121);
  VecDbl pt(2,0.0);
  LinearRegressionModelFactory lrmf;
  SurfpackModel* lrm = lrmf.Create(*sd);
  // First, see what the lrm is like with the default (no) scaling
  CPPUNIT_ASSERT(matches(0.0,(*lrm)(pt)));
  pt[0] = pt[1] = 2.0;
  CPPUNIT_ASSERT(matches(8.0,(*lrm)(pt)));
  StandardFitness sf;
  cout << "Standard Fitness: " << sf(*lrm,*sd) << endl;
  cout << lrm->asString() << endl; 
  cout << lrm->scaler()->asString() << endl; 
  delete sd;
  delete lrm;
  // Min is -2 and range is 4
}

