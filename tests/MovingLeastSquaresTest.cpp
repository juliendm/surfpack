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
#include "MovingLeastSquaresModel.h"
#include "MovingLeastSquaresTest.h"
#include "SurfpackModelTest.h"
#include "SurfData.h"
#include "SurfpackModel.h"
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

CPPUNIT_TEST_SUITE_REGISTRATION( MovingLeastSquaresTest );

void MovingLeastSquaresTest::setUp()
{

}

void MovingLeastSquaresTest::tearDown()
{

}

void MovingLeastSquaresTest::sineCurve()
{
  AxesBounds ab("-1 1 | -1 1");
  SurfData* sd = 0;
  VecUns gp(2,10);
  SurfpackInterface::CreateSample(sd,ab,gp);
  VecDbl responses(sd->size());
  for (unsigned i = 0; i < sd->size(); i++) {
    const VecDbl& pt = (*sd)(i);
    responses[i] = pt[0]*sin(20.0*pt[0])+pt[1]*sin(20.0*pt[1]);
  }
  sd->addResponse(responses);
  LRMBasisSet bs;
  bs.add("");
  bs.add("0");
  bs.add("1");
  bs.add("0 1");
  bs.add("0 0");
  bs.add("1 1");
  MovingLeastSquaresModel mlsm(*sd,bs);
  SurfpackModelTest::generalDerivativeTest(mlsm,ab);
   
}

void MovingLeastSquaresTest::generalModelTest()
{
  LRMBasisSet bs;
  bs.bases.push_back(VecUns());
  bs.bases.push_back(VecUns(2,0));
  bs.bases.push_back(VecUns(2,1));
  VecDbl mcfs(3,1.0);
  LinearRegressionModel my_model(2,bs,mcfs);
  AxesBounds ab("-2 2 | -2 2");
  SurfpackModelTest::generalDerivativeTest(my_model,ab);
}

