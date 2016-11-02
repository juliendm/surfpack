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
#include "ModelFactoryTest.h"
#include "SurfpackModel.h"
#include "MovingLeastSquaresModel.h"
#include "SurfpackInterface.h"
#include "SurfData.h"
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

CPPUNIT_TEST_SUITE_REGISTRATION( ModelFactoryTest );

void ModelFactoryTest::setUp()
{

}

void ModelFactoryTest::tearDown()
{

}
void ModelFactoryTest::simpleTest()
{
  SurfData* sd = SurfpackInterface::CreateSample("-2 2 | -2 2","10 10","sphere");
  //AxesBounds* ab = new AxesBounds("-2 2 | -2 2");
  //SurfData* sd = SurfpackInterface::CreateSample(ab, VecUns(10,10));
  SurfpackModelFactory* mlsf = new MovingLeastSquaresModelFactory;
  SurfpackModel* mlsm = mlsf->Create(*sd);
  VecDbl vd = surfpack::toVec<double>("0.0 0.0");
  cout << (*mlsm)(vd) << endl;
  delete mlsm;
  delete mlsf;
}

void ModelFactoryTest::argsTest()
{
  ParamMap args;
  args.insert(ModelParam(string("type"),string("polynomial")));
  args["order"] = "polynomial"; 
}

