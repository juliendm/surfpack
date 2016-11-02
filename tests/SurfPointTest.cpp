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

#include "SurfPointTest.h"
#include "SurfPoint.h"
#include "unittests.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::vector;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfPointTest );

void SurfPointTest::setUp()
{
  initialize();
  x1.resize(3);
  x1[0] = 3.0;
  x1[1] = -3.0;
  x1[2] = 0.0;

  x2.resize(1);
  x2[0] = 0.0;

  f1.resize(2);
  f1[0] = 1.0;
  f1[1] = -2.0;

  spPtr = new SurfPoint(x1);
  spPtr2 = new SurfPoint(x2, f1);
}

void SurfPointTest::tearDown()
{
  delete spPtr;
  delete spPtr2;
}

void SurfPointTest::testSurfPointPtrLessThan()
{
  SurfData::SurfPointSet s;
  s.insert(spPtr);
  s.insert(spPtr2);
  SurfData::SurfPointSet s2;
  s2.insert(spPtr2);
  s2.insert(spPtr);
}

void SurfPointTest::testConstructor()
{
  
  vector<double> pt;
  pt.push_back(1);
  SurfPoint s(pt);
  CPPUNIT_ASSERT_EQUAL(s.x[0], 1.0);
}

void SurfPointTest::testConstructorXSpecified()
{
  SurfPoint sp(x1);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], -3.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[2], 0.0);
  CPPUNIT_ASSERT(sp.x.size()==3);
  

}
void SurfPointTest::testConstructorXSpecifiedPlusOneF()
{
  SurfPoint sp(x2, 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 2.0);
  CPPUNIT_ASSERT(sp.x.size()==1);
  CPPUNIT_ASSERT(sp.f.size()==1);
}

void SurfPointTest::testConstructorXSpecifiedFVector()
{
  SurfPoint sp(x2, f1);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT(sp.x.size()==1);
  CPPUNIT_ASSERT(sp.f.size()==2);

}

void SurfPointTest::testConstructorFromIStreamBinary()
{
  // Create the binary SurfPoint in a file that will be read later
  //vector<double> xs(2);
  //xs[0] = 1.0;
  //xs[1] = 2.0;
  //vector<double> fs(2);
  //fs[0] = 3.0;
  //fs[1] = 4.0;
  //SurfPoint sp2(xs, fs);
  //ofstream outfile(fullPath("point1.sp"), ios::out | ios::binary);
  //sp2.writeBinary(outfile);
  //outfile.close();
 
  const string filename = "point1.sp";
  ifstream infile(filename.c_str());
  SurfPoint sp(2, 2, infile);  
  infile.close();
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 4.0);
  CPPUNIT_ASSERT(sp.x.size()==2);
  CPPUNIT_ASSERT(sp.f.size()==2);
}

void SurfPointTest::testConstructorFromIStreamText()
{
  const string filename = "point1.txt";
  ifstream infile(filename.c_str());
  string one_line;
  getline(infile,one_line);
  SurfPoint sp(2, 2, one_line, 0);  
  infile.close();
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 4.0);
  CPPUNIT_ASSERT(sp.x.size()==2);
  CPPUNIT_ASSERT(sp.f.size()==2);
}

void SurfPointTest::testCopyConstructor()
{
  SurfPoint sp2(x2, f1);
  SurfPoint sp(sp2);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT(sp.x.size()==1);
  CPPUNIT_ASSERT(sp.f.size()==2);
}

void SurfPointTest::testConstructorBadXSize()
{
  vector<double> x;
  SurfPoint sp(x); // should throw an exception
}

// Overloaded operators
void SurfPointTest::testOperatorAssignment()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(x1, 3.0);
  sp2 = sp;
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT(sp.x.size()==1);
  CPPUNIT_ASSERT(sp.f.size()==2);
}

void SurfPointTest::testOperatorAssignmentToSelf()
{
  SurfPoint sp(x2, f1);
  sp = sp;
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT(sp.x.size()==1);
  CPPUNIT_ASSERT(sp.f.size()==2);
}

void SurfPointTest::testOperatorEquality()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(sp);
  CPPUNIT_ASSERT(sp == sp2);
  CPPUNIT_ASSERT(sp.operator==(sp2));
  CPPUNIT_ASSERT(sp2.operator==(sp));
}

void SurfPointTest::testOperatorInequality()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(x1, 3.0);
  CPPUNIT_ASSERT(sp != sp2);
  CPPUNIT_ASSERT(sp.operator!=(sp2));
  CPPUNIT_ASSERT(sp2.operator!=(sp));
}

// Queries
void SurfPointTest::testXSize()
{
  CPPUNIT_ASSERT_EQUAL(spPtr->xSize(), static_cast<unsigned>(3));
  CPPUNIT_ASSERT_EQUAL(spPtr2->xSize(), static_cast<unsigned>(1));
}

void SurfPointTest::testFSize()
{
  CPPUNIT_ASSERT_EQUAL(spPtr->fSize(), static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(spPtr2->fSize(), static_cast<unsigned>(2));
}

void SurfPointTest::testX()
{
  vector<double> xvec = spPtr->X();
  CPPUNIT_ASSERT_EQUAL(xvec, x1);
}

void SurfPointTest::testFQuery()
{
  CPPUNIT_ASSERT_EQUAL(spPtr2->F(), 1.0);
  CPPUNIT_ASSERT_EQUAL(spPtr2->F(1), -2.0);
}

void SurfPointTest::testFQueryBadIndex()
{
  // should throw an exception, since spPtr has no response values
  double val = spPtr->F(6); 
}

// Commands
void SurfPointTest::testAddResponse()
{
  SurfPoint sp(*spPtr);
  unsigned newIndex = sp.addResponse();
  CPPUNIT_ASSERT_EQUAL(newIndex, static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(sp.fSize(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.F(), 0.0);
  newIndex = sp.addResponse(4.0);
  CPPUNIT_ASSERT_EQUAL(newIndex, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.fSize(), static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sp.F(1), 4.0);
}

void SurfPointTest::testFAssign()
{
  SurfPoint sp(*spPtr2);
  sp.F(1,3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 3.0);
}

void SurfPointTest::testFAssignBadIndex()
{
  SurfPoint sp(x1, f1);
  // Should throw an exception, since sp has no response values
  sp.F(3, 1.0); 
}

void SurfPointTest::testResize()
{
  SurfPoint x(vector<double>(1));
  x.resize(2);
  CPPUNIT_ASSERT(2==x.x.size());
  CPPUNIT_ASSERT_EQUAL(static_cast<unsigned>(2),x.xSize());
}

void SurfPointTest::testSetX()
{
  SurfPoint x4(*spPtr);
  x4.setX(2,1);
  CPPUNIT_ASSERT_EQUAL(1.0,x4.x[2]);  
  // Expand beyond current limits
  x4.setX(3,1);
  CPPUNIT_ASSERT_EQUAL(static_cast<unsigned>(4),x4.xSize());
  CPPUNIT_ASSERT(4==x4.x.size());
}
// I/O
void SurfPointTest::testWriteBinary()
{
  ofstream outfile(string("writePoint.sp").c_str(),ios::out | ios::binary);
  spPtr2->writeBinary(outfile);
  outfile.close();
  ifstream infile(string("writePoint.sp").c_str(),ios::in | ios::binary);
  SurfPoint sp(1, 2, infile);
  SurfPoint sp2(x2, f1);
  CPPUNIT_ASSERT(sp == sp2);
}

void SurfPointTest::testWriteText()
{
  ofstream outfile(string("writePoint.txt").c_str(),ios::out);
  spPtr2->writeText(outfile);
  outfile.close();
  ifstream infile(string("writePoint.txt").c_str(),ios::in);
  string one_line;
  getline(infile,one_line);
  SurfPoint sp(1, 2, one_line);
  SurfPoint sp2(x2, f1);
  CPPUNIT_ASSERT(sp == sp2);
}

void SurfPointTest::testReadBinary()
{
  const string filename = "point2.sp";
  ifstream infile(filename.c_str(), ios::in | ios::binary);
  SurfPoint sp(*spPtr2);
  sp.x[0] = 3.0;
  sp.readBinary(infile);
  CPPUNIT_ASSERT(sp == *spPtr2);
  infile.close();
}

void SurfPointTest::testReadText()
{
  const string filename = "point2.txt";
  ifstream infile(filename.c_str(), ios::in);
  SurfPoint sp(*spPtr2);
  sp.x[0] = 3.0;
  string one_line;
  getline(infile,one_line);
  sp.readText(one_line);
  CPPUNIT_ASSERT(sp == *spPtr2);
  infile.close();
}

void SurfPointTest::testStreamInsertion()
{
  // If this test doesn't throw an exception,
  // it is presumed to have worked
  ofstream blackhole("/dev/null",ios::out);
  blackhole << (*spPtr) << endl;
  cout << "End SurfPointTest" << endl;
}

