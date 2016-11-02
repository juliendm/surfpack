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
#include "KrigingModelTest.h"
#include "SurfpackModel.h"
#include "KrigingModel.h"
#include "DirectANNModel.h"
#include "unittests.h"
#include "SurfpackInterface.h"
#include "SurfData.h"
#include "surfpack.h"
#include "ModelFitness.h"
#include "surfpack.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::vector;
using std::ostream_iterator;
using std::ostringstream;

CPPUNIT_TEST_SUITE_REGISTRATION( KrigingModelTest );

void KrigingModelTest::setUp()
{

}

void KrigingModelTest::tearDown()
{

}
void KrigingModelTest::simpleTest()
{
  //SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","10 10","sphere");
  SurfData sd("data_samples.spd",3,1,0);
  cout << sd.size() << endl;
  KrigingModelFactory kmf;
  kmf.add("correlations","13.2792 3.07625 7.961");
  SurfpackModel* km = kmf.Build(sd);
  VecDbl responses = (*km)(sd);
  sd.addResponse(responses);
  sd.write("test_on_train.spd");
  
  //SurfData* sdp = SurfpackInterface::CreateSample("2.53e06 2.55e06 | 7.64e05 7.66e05 | 148 149","10 10 10","");
  SurfData* sdp = SurfpackInterface::CreateSample("1.53e06 2.55e06 | 7.04e05 8.66e05 | 140 159","10 10 10","");
  VecDbl responses2 = (*km)(*sdp);
  sdp->addResponse(responses2);
  sdp->write("test_data.spd");
  delete sdp;
}

