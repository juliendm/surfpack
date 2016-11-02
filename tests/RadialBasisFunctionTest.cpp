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
#include "SurfpackMatrix.h"
#include "surfpack.h"
#include "RBFNetSurface.h"
#include "RadialBasisFunctionModel.h"
#include "RadialBasisFunctionTest.h"
#include "SurfData.h"
#include "SurfpackModel.h"
#include "SurfpackModelTest.h"
#include "SurfpackInterface.h"
#include "KrigingModel.h"
#include "DirectANNModel.h"
#include "AxesBounds.h"
#include "unittests.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::vector;
using std::ostream_iterator;
using std::ostringstream;
using surfpack::shared_rng;

CPPUNIT_TEST_SUITE_REGISTRATION( RadialBasisFunctionTest );

void RadialBasisFunctionTest::setUp()
{

}

void RadialBasisFunctionTest::tearDown()
{

}

void RadialBasisFunctionTest::generalModelTest()
{
  VecRbf rbfs;
  rbfs.push_back(RadialBasisFunction("0.0 0.0","1.0 1.0"));
  VecDbl cfs(1,1.0);
  RadialBasisFunctionModel rbf_model(rbfs,cfs);
  AxesBounds ab("-2 2 | -2 2");
  cout << rbf_model.asString() << endl;
  //SurfpackModelTest::generalDerivativeTest(rbf_model,ab);
}
void RadialBasisFunctionTest::hillValleyTest()
{
  VecRbf rbfs;
  //rbfs.push_back(RadialBasisFunction("0.5 0.5","0.5 0.5"));
  //rbfs.push_back(RadialBasisFunction("-0.5 -0.5","0.5 0.5"));
  rbfs.push_back(RadialBasisFunction("0.5 0.5","10 1"));
  //rbfs.push_back(RadialBasisFunction("-0.5 -0.5","1 1"));
  VecDbl cfs(1,1.0);
  //cfs[1] = -1.0;
  RadialBasisFunctionModel rbf_model(rbfs,cfs);
  AxesBounds ab("-2 2 | -2 2");
  cout << rbf_model.asString() << endl;
  SurfpackModelTest::generalDerivativeTest(rbf_model,ab);
}

void RadialBasisFunctionTest::eTest()
{
  return;
  VecRbf rbfs;
  VecDbl center(2);
  VecDbl radius(2);
  VecDbl cfs;
  for (unsigned i = 0; i < 100; i++) {
    center[0] = shared_rng().rand()*20.0-10.0;
    center[1] = shared_rng().rand()*20.0-10.0;
    radius[0] = shared_rng().rand()*10.0;
    radius[1] = shared_rng().rand()*10.0;
    rbfs.push_back(RadialBasisFunction(center,radius));
    cfs.push_back(shared_rng().rand()*5.0-2.5);
  }
  RadialBasisFunctionModel rbf_model(rbfs,cfs);
  AxesBounds ab("-30 30 | -30 30");
  cout << rbf_model.asString() << endl;
  VecUns grid_points(2,500);
  SurfData* sd = 0;
  SurfpackInterface::CreateSample(sd,ab,grid_points);
  SurfpackModel& sm = rbf_model;
  rbf_model(center); 
  time_t elapsed = -time(0);
  VecDbl v = sm(*sd);
  elapsed += time(0);
  cout << "elapsed: " << elapsed << endl;
}

void RadialBasisFunctionTest::partitionTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2 | -2 2","10 10 10","sphere");
  RBFNetSurface* rbf = new RBFNetSurface(sd); 
  rbf->partition(*sd,25);
  delete sd; sd = 0;
}

void RadialBasisFunctionTest::centroidTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","10 10","sphere");
  SurfPoint sp = computeCentroid(*sd); 
  delete sd; sd = 0;
}

void RadialBasisFunctionTest::centroidTest2()
{
  SurfData sd;
  sd.addPoint(SurfPoint(surfpack::toVec<double>("1 2 3")));
  sd.addPoint(SurfPoint(surfpack::toVec<double>("4 2 9")));
  sd.addPoint(SurfPoint(surfpack::toVec<double>("-2 2 3")));
  SurfPoint c = computeCentroid(sd);
  CPPUNIT_ASSERT(matches(c[0],1.0));
  CPPUNIT_ASSERT(matches(c[1],2.0));
  CPPUNIT_ASSERT(matches(c[2],5.0));
}

void RadialBasisFunctionTest::updateCentroidTest()
{
  SurfData sd;
  sd.addPoint(SurfPoint(surfpack::toVec<double>("1 2 3")));
  sd.addPoint(SurfPoint(surfpack::toVec<double>("4 2 9")));
  sd.addPoint(SurfPoint(surfpack::toVec<double>("-2 2 3")));
  VecDbl centroid(3,0.0);
  updateCentroid(centroid,sd(0),0);
  updateCentroid(centroid,sd(1),1);
  updateCentroid(centroid,sd(2),2);
  VecDbl c = centroid;
  CPPUNIT_ASSERT(matches(c[0],1.0));
  CPPUNIT_ASSERT(matches(c[1],2.0));
  CPPUNIT_ASSERT(matches(c[2],5.0));

}

void RadialBasisFunctionTest::cvtTest()
{
  AxesBounds ab("-2 2 | -2 2");
  SurfData sd = cvts(ab);
  SurfData radiuses = radii(sd);
  sd.write("cvts.spd");
}

void RadialBasisFunctionTest::createTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","10 10","moderatepoly");
  RadialBasisFunctionModelFactory rbfmf;
  SurfpackModel* model = rbfmf.Create(*sd);
  VecDbl est = (*model)(*sd);
  sd->addResponse(est,"est");
  sd->write("rbftest.spd");
  delete sd;
  delete model;
}

