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

#include "SurfpackCommonTest.h"
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"
#include "SurfpackInterpreter.h"
#include "unittests.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::vector;
using std::ostream_iterator;


CPPUNIT_TEST_SUITE_REGISTRATION( SurfpackCommonTest );

void SurfpackCommonTest::setUp()
{
}

void SurfpackCommonTest::tearDown()
{
}

void SurfpackCommonTest::testConstructor()
{
  CPPUNIT_ASSERT( true );
}


void SurfpackCommonTest::matrixMultiplyTest()
{
  SurfpackMatrix< double > m1(2,3,"1 0 2 1 1 0",true);
  SurfpackMatrix< double > result;
  surfpack::matrixMatrixMult(result,m1,m1,'N','T');
  CPPUNIT_ASSERT(matches(result(0,0),6.0));
  CPPUNIT_ASSERT(matches(result(0,1),2.0));
  CPPUNIT_ASSERT(matches(result(1,0),2.0));
  CPPUNIT_ASSERT(matches(result(1,1),1.0));
}

void SurfpackCommonTest::matrixVectorMultTest()
{

  SurfpackMatrix< double > m1(2,3,"1 0 2 1 1 0",true);
  vector< double > c(3);
  c[0] = 1;
  c[1] = 2;
  c[2] = 3;
  vector< double > result;
  surfpack::matrixVectorMult(result,m1,c);
  CPPUNIT_ASSERT(matches(result[0],8.0));
  CPPUNIT_ASSERT(matches(result[1],2.0));
  c.resize(2);
  c[0] = 1;
  c[1] = 2;
  surfpack::matrixVectorMult(result,m1,c,'T');
  CPPUNIT_ASSERT(matches(result[0],1.0));
  CPPUNIT_ASSERT(matches(result[1],4.0));
  CPPUNIT_ASSERT(matches(result[2],1.0));
  copy(result.begin(),result.end(),ostream_iterator<double>(cout,"\n"));
}

void SurfpackCommonTest::matrixVectorMultTransTest()
{

  SurfpackMatrix< double > m1(3,2,"1 0 2 1 1 0",true);
  vector< double > c(3);
  c[0] = 1;
  c[1] = 2;
  c[2] = 3;
  vector< double > result;
  surfpack::matrixVectorMult(result,m1,c,'T');
  copy(result.begin(),result.end(),ostream_iterator<double>(cout,"\n"));
  CPPUNIT_ASSERT(matches(result[0],7.0));
  CPPUNIT_ASSERT(matches(result[1],3.0));
}

void SurfpackCommonTest::matrixReshapeTest()
{
  SurfpackMatrix<double> sm(2, 3, "1 2 3 4 5 6",true);
  CPPUNIT_ASSERT(matches(sm(0,0),1.0));
  CPPUNIT_ASSERT(matches(sm(1,0),2.0));
  cout << sm.asString() << endl;
  sm.resize(4,5);
  CPPUNIT_ASSERT(matches(sm(0,0),1.0));
  CPPUNIT_ASSERT(matches(sm(1,2),6.0));
  CPPUNIT_ASSERT(matches(sm(0,1),3.0));
  CPPUNIT_ASSERT(matches(sm(2,0),0.0));
  CPPUNIT_ASSERT(matches(sm(3,0),0.0));
  cout << sm.asString() << endl;
}

void SurfpackCommonTest::matrixResizeCTest()
{
  SurfpackMatrix<double> sm(2, 3, "1 2 3 4 5 6",false);
  CPPUNIT_ASSERT(matches(sm(0,0),1.0));
  cout << sm.asString() << endl;
  sm.resize(4,5);
  CPPUNIT_ASSERT(matches(sm(0,0),1.0));
  cout << sm.asString() << endl;

}

void SurfpackCommonTest::stringToVecUnsTest()
{
  string s("1 2 3");
  vector<unsigned> v = surfpack::toVec<unsigned>(s);
  copy(v.begin(),v.end(),ostream_iterator<unsigned>(cout,"\n"));
  CPPUNIT_ASSERT(v.size() == 3);
  CPPUNIT_ASSERT(v[0] == 1);
  CPPUNIT_ASSERT(v[1] == 2);
  CPPUNIT_ASSERT(v[2] == 3);
  string s2("");
  vector<unsigned> v2 = surfpack::toVec<unsigned>(s2);
  CPPUNIT_ASSERT(v2.empty());
  
}

void SurfpackCommonTest::weightedAvgTest()
{
  VecDbl x(3,1.0);
  VecDbl y(3,0.0);
  VecDbl z = surfpack::weightedAvg(x,y,.4);
  CPPUNIT_ASSERT(matches(z[0],.4));
}

void SurfpackCommonTest::toString()
{
}

void SurfpackCommonTest::fromVec()
{

  std::vector<unsigned> uv(4,3);
  CPPUNIT_ASSERT(surfpack::fromVec<unsigned>(uv) == string("3 3 3 3"));
  
  //cout << surfpack::fromVec<unsigned>(uv) << endl;
  std::vector<std::string> sv;
  sv.push_back("sphere");
  sv.push_back("rosenbrock");
  sv.push_back("rastrigin");
  CPPUNIT_ASSERT(surfpack::fromVec<std::string>(sv) == string("sphere rosenbrock rastrigin"));
  //CPPUNIT_ASSERT(surfpack::fromVec<std::string>(sv) == string("sphere rosenbrock rastrigin"));

}

void SurfpackCommonTest::blockTests()
{
  unsigned p = 5;
  unsigned n = 17;
  CPPUNIT_ASSERT(surfpack::block_low(0,p,n) == 0);
  CPPUNIT_ASSERT(surfpack::block_low(1,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_low(2,p,n) == 6);
  CPPUNIT_ASSERT(surfpack::block_low(3,p,n) == 10);
  CPPUNIT_ASSERT(surfpack::block_low(4,p,n) == 13);
  CPPUNIT_ASSERT(surfpack::block_low(0,p,n) == 0);
  CPPUNIT_ASSERT(surfpack::block_high(0,p,n) == 2);
  CPPUNIT_ASSERT(surfpack::block_high(1,p,n) == 5);
  CPPUNIT_ASSERT(surfpack::block_high(2,p,n) == 9);
  CPPUNIT_ASSERT(surfpack::block_high(3,p,n) == 12);
  CPPUNIT_ASSERT(surfpack::block_high(4,p,n) == 16);
  CPPUNIT_ASSERT(surfpack::block_size(0,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_size(1,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_size(2,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_size(3,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_size(4,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_owner(0,p,n) == 0);
  CPPUNIT_ASSERT(surfpack::block_owner(1,p,n) == 0);
  CPPUNIT_ASSERT(surfpack::block_owner(2,p,n) == 0);
  CPPUNIT_ASSERT(surfpack::block_owner(3,p,n) == 1);
  CPPUNIT_ASSERT(surfpack::block_owner(4,p,n) == 1);
  CPPUNIT_ASSERT(surfpack::block_owner(5,p,n) == 1);
  CPPUNIT_ASSERT(surfpack::block_owner(6,p,n) == 2);
  CPPUNIT_ASSERT(surfpack::block_owner(7,p,n) == 2);
  CPPUNIT_ASSERT(surfpack::block_owner(8,p,n) == 2);
  CPPUNIT_ASSERT(surfpack::block_owner(9,p,n) == 2);
  CPPUNIT_ASSERT(surfpack::block_owner(10,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_owner(11,p,n) == 3);
  CPPUNIT_ASSERT(surfpack::block_owner(12,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_owner(13,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_owner(14,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_owner(15,p,n) == 4);
  CPPUNIT_ASSERT(surfpack::block_owner(16,p,n) == 4);
}

