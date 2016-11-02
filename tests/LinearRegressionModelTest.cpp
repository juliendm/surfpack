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

#include "LinearRegressionModelTest.h"
#include "LinearRegressionModel.h"
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
using surfpack::shared_rng;

const int TESTDIMS = 2;
const int GRIDSIZE = 50;
CPPUNIT_TEST_SUITE_REGISTRATION( LinearRegressionModelTest );


// BMA: These functions previously in ModelFitness.h, but not used there
template< typename T>
std::vector< T >& vecSubInPlace(std::vector< T >& sub, std::vector< T >& min)
{
  assert(sub.size() == min.size());
  typedef typename std::vector< T >::iterator VecIt;
  VecIt subIt, minIt;
  for (subIt = sub.begin(), minIt = min.begin();
	subIt != sub.end();
	++subIt, ++minIt) {
    *subIt -= *minIt;
  }
  return sub;
}

template< typename T>
std::vector< T > vecSub(std::vector< T >& sub, std::vector< T >& min)
{
  assert(sub.size() == min.size());
  std::vector< T > result(sub.size());
  typedef typename std::vector< T >::iterator VecIt;
  VecIt subIt, minIt, resIt;
  for (subIt = sub.begin(), minIt = min.begin(), resIt = result.begin();
	subIt != sub.end();
	++subIt, ++minIt, ++resIt) {
    *resIt = *subIt - *minIt;
  }
  return result;
}


void generateDerivPlot(const std::string& plotname, const std::string& datafilename, int x, int y1, int y2)
{

  std::ofstream outfile("gnuplotscript",std::ios::out);
  outfile << "set term png\n"
	     "set output \"" << plotname << ".png\"\n"
	     "set autoscale\n"
             "plot  \"" << datafilename << ".spd\" using "
	  << x << ":"  << y1
	 << ", \"" << datafilename << ".spd\" using " 
	  << x << ":"  << y2
	 << ", \"" << datafilename << "pt\" using 1:2 points 10\n" ;
  outfile.close();
  system("gnuplot ./gnuplotscript\n");
  system("rm gnuplotscript");

}

void LinearRegressionModelTest::setUp()
{
  // Create data set
  ab = new AxesBounds(string("-2 2 | -2 2"));
  vector< unsigned > grid_points(2);
  grid_points[0] = grid_points[1] = GRIDSIZE;
  sd = 0;
  randsd = 0;
  vector< string > test_functions;
  // Points on grid
  SurfpackInterface::CreateSample(sd,*ab,grid_points,test_functions);
  // Random sample of points
  SurfpackInterface::CreateSample(randsd,*ab,10,test_functions);

}

void LinearRegressionModelTest::tearDown()
{
  delete sd;
  delete ab;
  delete randsd;
}

void LinearRegressionModelTest::constructorTest()
{
  VecDbl cf(1,1.0);
  LRMBasisSet bs;
  bs.bases.push_back(VecUns());
  LinearRegressionModel lrm(1,bs,cf);
  CPPUNIT_ASSERT(matches(cf[0],lrm.coeffs[0]));
  CPPUNIT_ASSERT(lrm.size() == 1);

}

void LinearRegressionModelTest::unityBasisTest()
{
  VecUns list;
  bs.bases.push_back(list);
  VecDbl x(1,0.0);
  VecUns vars(1,0);
  CPPUNIT_ASSERT(matches(bs.eval(0,x),1.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),0.0));
}

void LinearRegressionModelTest::singleLinearTest()
{

  VecUns list(1,0);
  bs.bases.push_back(list);
  VecDbl x(1,0.0);
  VecUns vars(1,0);
  // The term is x0
  // At x0 = 0, eval = 0 and deriv = 1
  CPPUNIT_ASSERT(matches(bs.eval(0,x),0.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),1.0));
  x[0] = 1.0;
  // At x0 = 1, eval = 1 and deriv = 1
  CPPUNIT_ASSERT(matches(bs.eval(0,x),1.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),1.0));
  x[0] = 2.0;
  // At x0 = 2, eval = 2 and deriv = 1
  CPPUNIT_ASSERT(matches(bs.eval(0,x),2.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),1.0));
  // Now differentiate again with respect to x0
  vars.push_back(0);
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),0.0));
}

void LinearRegressionModelTest::singleQuadraticTest()
{
  VecUns list(2,0);
  bs.bases.push_back(list);
  VecDbl x(1,0.0);
  VecUns vars(1,0);
  // The term is x0^2
  // At x0 = 0, eval = 0 and deriv = 0
  CPPUNIT_ASSERT(matches(bs.eval(0,x),0.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),0.0));
  x[0] = 1.0;
  // At x0 = 1, eval = 1 and deriv = 2
  CPPUNIT_ASSERT(matches(bs.eval(0,x),1.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),2.0));
  x[0] = 2.0;
  // At x0 = 2, eval = 4 and deriv = 4
  CPPUNIT_ASSERT(matches(bs.eval(0,x),4.0));
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),4.0));
  // Now differentiate again with respect to x0
  vars.push_back(0);
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),2.0));
  // Now differentiate again with respect to x0
  vars.push_back(0);
  CPPUNIT_ASSERT(matches(bs.deriv(0,x,vars),0.0));
}

void LinearRegressionModelTest::lineEvalTest()
{
  VecUns list;
  bs.bases.push_back(list);
  list.push_back(0);
  bs.bases.push_back(list);
  VecDbl cf(2,1.0);
  LinearRegressionModel lrm(1,bs,cf); // lrm = x[0] + 1;
  VecDbl x(1,1.0);
  CPPUNIT_ASSERT(lrm.size() == 1);
  CPPUNIT_ASSERT(matches(lrm(x),2.0));
  VecDbl g = lrm.gradient(x);
  CPPUNIT_ASSERT(matches(g[0],1.0));

  // Evaluate at x[0] = -3.0;
  x[0] = -3.0;
  CPPUNIT_ASSERT(matches(lrm(x),-2.0));
  VecDbl g1 = lrm.gradient(x);
  CPPUNIT_ASSERT(matches(g1[0],1.0));

  // Change function to 2*x[0] + 3;
  lrm.coeffs[0] = 3.0;
  lrm.coeffs[1] = 2.0;
  // Evaluation = 2*-3 + 3;
  CPPUNIT_ASSERT(matches(lrm(x),-3.0));
  VecDbl g2 = lrm.gradient(x);
  cout << "grad: " << g2[0] << endl;
  CPPUNIT_ASSERT(matches(g2[0],2.0));
}

void LinearRegressionModelTest::quadratic2DTest()
{
  VecUns list;
  bs.bases.push_back(list); // Unity
  list.push_back(0);
  bs.bases.push_back(list); // x[0]
  list.push_back(1);
  bs.bases.push_back(list); // x[0]*x[1]
  list[0] = 1;
  bs.bases.push_back(list); // *x[1]*x[1]

  VecDbl cf(4,-1.0);
  cf[1] = 2.0;
  cf[2] = -3.0;
  cf[3] = 4.0;
  // lrm = -1 + 2.0*x[0] - 3.0*x[0]*x[1] + 4.0*x[1]^2;
  LinearRegressionModel lrm(2,bs,cf); 
  VecDbl x(2,1.0);
  // lrm(1,1) = -1+2(1)-3(1)(1)+4(1)(1) = 2
  CPPUNIT_ASSERT(lrm.size() == 2);
  CPPUNIT_ASSERT(matches(lrm(x),2.0));
  VecDbl g = lrm.gradient(x);
  // dlrm/dx0 = 2.0 - 3.0*x[1]
  // dlrm/dx1 = -3.0*x[0]+8.0*x[1]
  CPPUNIT_ASSERT(matches(g[0],-1.0));
  CPPUNIT_ASSERT(matches(g[1],5.0));

   //Evaluate at (-3,-1)
  x[0] = -3.0;
  x[1] = -1.0;
  // -1+2(-3)-3(-3)(-1)+4(-1)(-1) = -12
  CPPUNIT_ASSERT(matches(lrm(x),-12.0));
  VecDbl g1 = lrm.gradient(x);
  CPPUNIT_ASSERT(matches(g1[0],5.0));
  CPPUNIT_ASSERT(matches(g1[1],1.0));
}

void LinearRegressionModelTest::plotTest1()
{
  

  VecUns list;
  bs.bases.push_back(list); // Unity
  list.push_back(0);
  bs.bases.push_back(list); // x[0]
  list.push_back(1);
  bs.bases.push_back(list); // x[0]*x[1]
  list[0] = 1;
  bs.bases.push_back(list); // *x[1]*x[1]

  VecDbl cf(4,-1.0);
  cf[1] = 2.0;
  cf[2] = -3.0;
  cf[3] = 4.0;
  // lrm = -1 + 2.0*x[0] - 3.0*x[0]*x[1] + 4.0*x[1]^2;
  LinearRegressionModel lrm(2,bs,cf); 

  // Prepare a basis set for a hyperplane model 
  LRMBasisSet hpbs; // hyperplane basis set
  list.resize(1);
  hpbs.bases.push_back(VecUns()); // Unity
  for (unsigned i = 0; i < randsd->xSize(); i++) {
    list[0] = i;
    hpbs.bases.push_back(list);
  }


  SurfData& rsd = *randsd;
  for (unsigned i = 0; i < rsd.size(); i++) {
    VecDbl pt = rsd(i);
    double resp = lrm(pt);
    VecDbl grad = lrm.gradient(pt);
    // solve for b in y = m1x1 + m2x2 + ... + b;
    double b = resp;
    for (unsigned j = 0; j < pt.size(); j++) {
      b -= pt[j]*grad[j];
    }

    // Now create a tangent hyperplane to the surface
    VecDbl hyp_coeff = grad;
    reverse(hyp_coeff.begin(),hyp_coeff.end());
    hyp_coeff.push_back(b);
    reverse(hyp_coeff.begin(),hyp_coeff.end());
    LinearRegressionModel tang_plane(pt.size(),hpbs,hyp_coeff);

    // Now prepare a set of axes, so that we can plot a slice of both
    // functions along each dimension
    // This is a set of axes along which each dimension is fixed
    vector<AxesBounds::Axis> axes(pt.size());
    for (unsigned j = 0; j < pt.size(); j++) {
      axes[j].min = axes[j].max = pt[j];
    }
    VecUns dims(pt.size(),1);
    
    for (unsigned dim = 0; dim < pt.size(); dim++) {
      vector<AxesBounds::Axis> var_axes = axes;
      var_axes[dim].min = -2.0;
      var_axes[dim].max = 2.0;
      AxesBounds var_ab(var_axes);
      VecUns var_dims = dims;
      var_dims[dim] = GRIDSIZE;
      SurfData* test_data = 0;
      vector<string> test_functions;
      SurfpackInterface::CreateSample(test_data,var_ab,var_dims,test_functions);
      SurfData& td = *test_data;
      // Now evaluate both models on this data
      ///\\todo Creat a base model class that can do operator()(SurfData&)
      VecDbl lrm_resp(td.size());
      VecDbl hyp_resp(td.size());
      for (unsigned k = 0; k < td.size(); k++) {
        lrm_resp[k] = lrm(td(k));
        hyp_resp[k] = tang_plane(td(k));
      }
      unsigned y1 = td.addResponse(lrm_resp) + td.xSize() + 1; // gnuplot counts cols from 1
      unsigned y2 = td.addResponse(hyp_resp) + td.xSize() + 1;
      unsigned xdim = dim + 1;
      string s1("deriv");
      ostringstream os;
      os << s1 << i << "_" << dim ;
      string plotfilename = os.str();
      string datafilename = plotfilename + ".spd";
      string ptfilename = plotfilename + "pt";
      std::ofstream ptfile(ptfilename.c_str(),ios::out);
      ptfile << pt[dim] << " " << resp << endl;
      ptfile.close(); 
      td.write(datafilename);
      generateDerivPlot(plotfilename,plotfilename,xdim,y1,y2); 
      cout << "Pt: " << i << " dim: " << dim << endl;
      cout << lrm.asString() << tang_plane.asString() << endl;
      delete test_data;
    }
    
  } 
}

void LinearRegressionModelTest::termPrinterTest()
{
  LinearRegressionModelFactory::CreateLRM(2,2);
  VecUns u, v;
  u.push_back(1);
  u.push_back(10);
  u.push_back(100);
  v.push_back(0);
  v.push_back(5);
  v.push_back(30);
  vecSub<unsigned>(u,v);
  copy(u.begin(),u.end(),std::ostream_iterator<unsigned>(std::cout,"\n"));  
}

void LinearRegressionModelTest::createModelTest()
{
  SurfData* sd = SurfpackInterface::CreateSample(string("-2 2 | -2 2"),
    string("8 8"),string("sphere"));
  LinearRegressionModelFactory lrmf;
  SurfpackModel* lrm = lrmf.Create(*sd);
  cout << lrm->asString() << endl;
  return;
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[0],0.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[1],0.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[2],1.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[3],0.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[4],0.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<LinearRegressionModel*>(lrm)->coeffs[5],1.0));
  delete lrm;
  delete sd;
}

//extern "C" double gsl_ran_fdist_pdf(double,double,double);
//void LinearRegressionModelTest::FTest()
//{
//  // Make a set of data that is a full quadratic fit plus some random noise
//  SurfData* sd = SurfpackInterface::CreateSample(string("-2 2 | -2 2"),
//    string("8 8"),string("sphere"));
//  LRMBasisSet bs = LinearRegressionModelFactory::CreateLRM(2,2);
//  CPPUNIT_ASSERT(bs.size() == 6);
//  VecDbl coeffs = surfpack::toVec<double>(string("1 -2 3 -1 3 1"));
//  LinearRegressionModel lrmTemp(2,bs,coeffs);
//  VecDbl responses = lrmTemp(*sd);
//  // Add some random noise to the data
//  for (unsigned i = 0; i < responses.size(); i++) {
//    responses[i] += (shared_rng().rand()-0.5)/2.0;
//  }
//  unsigned new_index = sd->addResponse(responses);
//  sd->setDefaultIndex(new_index);
//  // Now fit a full model to this data
//  LinearRegressionModelFactory lrmf;
//  LinearRegressionModel* lrmFull = lrmf.Create(*sd);
//  cout << lrmFull->asString() << endl;
//  StandardFitness sf;
//  double ssr_full = sf(*lrmFull,*sd);
//  // Now remove a term to create a reduced model
//  bs.bases.erase(bs.bases.begin()+2);
//  VecDbl coeffs_reduced = LinearRegressionModel::lrmSolve(bs,ScaledSurfData(NonScaler(),*sd));
//  LinearRegressionModel lrmReduced(2,bs,coeffs_reduced);
//  double ssr_reduced = sf(lrmReduced,*sd);
//  cout << "Fitness Full: " << ssr_full << endl;
//  cout << "Fitness Reduced: " << ssr_reduced << endl;
//
//  // compute degrees of freedom
//  unsigned df_num = 1;
//  unsigned df_denom = sd->size() - coeffs.size();
//  double Fstat = (ssr_reduced - ssr_full)*df_denom/ssr_full;
//  double fpdf = gsl_ran_fdist_pdf(Fstat,df_num,df_denom);
//  cout << "Fstat: " << Fstat << " fpdf: " << fpdf << endl;
//  delete lrmFull;
//}

