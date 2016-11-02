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
#include "SurfpackModelTest.h"
#include "SurfpackModel.h"
#include "KrigingModel.h"
#include "DirectANNModel.h"
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"
#include "SurfpackInterface.h"
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

const int TESTDIMS = 2;
const int GRIDSIZE = 50;
CPPUNIT_TEST_SUITE_REGISTRATION( SurfpackModelTest );

void generateDerivPlot1(const std::string& plotname, const std::string& datafilename, int x, int y1, int y2)
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

void SurfpackModelTest::setUp()
{
  // Create data set
  //ab = new AxesBounds(string("-2 2 | -2 2"));
  ab = new AxesBounds(string("-1 1 | -1 1"));
  vector< unsigned > grid_points(2, GRIDSIZE);
  sd = 0;
  randsd = 0;
  vector< string > test_functions;
  // Points on grid
  SurfpackInterface::CreateSample(sd,*ab,grid_points,test_functions);
  // Random sample of points
  SurfpackInterface::CreateSample(randsd,*ab,10,test_functions);

}

void SurfpackModelTest::tearDown()
{
  delete sd;
  delete ab;
  delete randsd;
}


void SurfpackModelTest::plotTest1()
{
  
  
// Kriging Specific stuff
  VecVecDbl bases;
  VecDbl corr(2);
  corr[0] = 1.0;
  corr[1] = 10.0;

  VecDbl pt(2);
  pt[0] = 0.0;
  pt[1] = 0.0;
  bases.push_back(pt);
  pt[0] = 1.0;
  pt[1] = 1.0;
  bases.push_back(pt);
  pt[0] = -1.0;
  pt[1] = -1.0;
  bases.push_back(pt);
  pt[0] = -1.5;
  pt[1] = -1.5;
  bases.push_back(pt);
  KrigingBasisSet bs(bases,corr);

  VecDbl cf(4,-1.0);
  cf[1] = 2.0;
  cf[2] = -3.0;
  cf[3] = 4.0;
  // lrm = -1 + 2.0*x[0] - 3.0*x[0]*x[1] + 4.0*x[1]^2;
  //KrigingModel lrm(bs,cf); 

// ANN Direct stuff
  const int neurons = 4;
  const int ndims = 2;
  MtxDbl weights(neurons,ndims+1,true);
  VecDbl oweights(neurons+1);
  for (unsigned i = 0; i < weights.getNRows(); i++) {
    oweights[i] = shared_rng().rand();
    for (unsigned j = 0; j < weights.getNCols(); j++) {
      weights(i,j) = shared_rng().rand();
    }
  }
  oweights.back() = shared_rng().rand();
  DirectANNBasisSet bsweights(weights);
  DirectANNModel lrm(bsweights,oweights);

  // Prepare a basis set for a hyperplane model 
  LRMBasisSet hpbs; // hyperplane basis set
  VecUns list(1,0);
  hpbs.bases.push_back(VecUns()); // Unity
  for (unsigned i = 0; i < randsd->xSize(); i++) {
    list[0] = i;
    hpbs.bases.push_back(list);
  }


  SurfData& rsd = *randsd;
  for (unsigned i = 0; i < rsd.size(); i++) {
    VecDbl pt = rsd(i);
    double resp = lrm(pt);
    cout << "Response: " << resp << endl;
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
      var_axes[dim].minIsMax = false;
      AxesBounds var_ab(var_axes);
      VecUns var_dims = dims;
      var_dims[dim] = GRIDSIZE;
      vector<string> test_functions;
      SurfData* test_data = var_ab.sampleGrid(var_dims,test_functions);
      //cout << "var dims: ";
      //copy(var_dims.begin(),var_dims.end(),std::ostream_iterator<double>(cout," "));
      //cout << "\nVar ab\n" << var_ab.asString() << endl;
      //continue;
      //SurfpackInterface::CreateSample(test_data,var_ab,var_dims,test_functions);
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
      generateDerivPlot1(plotfilename,plotfilename,xdim,y1,y2); 
      cout << "Pt: " << i << " dim: " << dim << endl;
      cout << lrm.asString() << tang_plane.asString() << endl;
      delete test_data;
    }
    
  } 
}

void SurfpackModelTest::graphicalDerivTest()
{

}

void SurfpackModelTest::modelSample(const SurfpackModel& model)
{
 
  VecDbl resps = model(*sd);
  copy(resps.begin(),resps.end(),std::ostream_iterator<double>(cout,"\n"));
  double minimum1 = *(min_element(resps.begin(),resps.end()));
  double maximum1 = *(max_element(resps.begin(),resps.end()));
  printf("min: %f max: %f\n",minimum1,maximum1);
}

void SurfpackModelTest::modelSampleTest()
{
  
  // Create model;
  const int neurons = 4;
  const int ndims = 2;
  MtxDbl weights(neurons,ndims+1,true);
  VecDbl oweights(neurons+1);
  for (unsigned i = 0; i < weights.getNRows(); i++) {
    oweights[i] = (shared_rng().rand())*.1-.05;
    for (unsigned j = 0; j < weights.getNCols(); j++) {
      weights(i,j) = (shared_rng().rand())*.1-.05;
    }
  }
  oweights.back() = (shared_rng().rand())*.1-.05;
  DirectANNBasisSet bsweights(weights);
  DirectANNModel lrm(bsweights,oweights);

  
  modelSample(lrm);
  return;


  // Prepare a basis set for a hyperplane model 
  LRMBasisSet hpbs; // hyperplane basis set
  VecUns list(1,0);
  hpbs.bases.push_back(VecUns()); // Unity
  for (unsigned i = 0; i < randsd->xSize(); i++) {
    list[0] = i;
    hpbs.bases.push_back(list);
  }
  VecDbl cofs(randsd->xSize()+1,1.0);
  LinearRegressionModel lrm2(randsd->xSize(),hpbs,cofs);
  modelSample(lrm2);
}

  const double w0 = .3;
  const double w1 = -.2;
  const double w2 = -.02;
  const double w00 = .15;
  const double w01 = .1;
  const double w02 = -.05;
  const double w10 = .09;
  const double w11 = -.02;
  const double w12 = .03;
double func1(double x0, double x1,bool pderiv)
{
  double y,s,g,h;
  g = w00*x0+w01*x1+w02;
  h = w10*x0+w11*x1+w12;
  s = w0*tanh(g) + w1*tanh(h) + w2;
  y = tanh(s);
  //printf("g %f h %f s %f y%f\n",g,h,s,y);
  printf("y(%f,%f): %f\n",x0,x1,y);

  if (pderiv) {
    double dyd0, dyds, dsd0, dsdg, dsdh, dgd0, dhd0;
    dhd0 = w10;
    dgd0 = w00;
    dsdg = w0*(1-tanh(g)*tanh(g));
    dsdh = w1*(1-tanh(h)*tanh(h));
    dyds = 1-tanh(s)*tanh(s);
    dsd0 = dsdg*dgd0+dsdh*dhd0;
    dyd0 = dyds*dsd0;
    printf("dy(.5,.2)/d0: %f\n",dyd0);
  }
  return y;
}

void SurfpackModelTest::manualANNTest()
{
  double x0 = .5;
  double x1 = .4;
  double y = func1(x0,x1,true);
  
  double epsilon = 1.0;
  do {
    double xp = x0+epsilon;
    double xm = x0-epsilon;
    double dy = func1(xp,x1,false)-func1(xm,x1,false);
    double dx0 = xp-xm;
    double numder = dy/dx0;
    printf("eps: %f slope: %f\n",epsilon,numder);
    epsilon *= .5;
  } while (epsilon > 1e-3);

  // Now see what the model says
  MtxDbl wghts(2,3,true);
  wghts(0,0) = w00;
  wghts(0,1) = w01;
  wghts(0,2) = w02;
  wghts(1,0) = w10;
  wghts(1,1) = w11;
  wghts(1,2) = w12;
  VecDbl ow(3);
  ow[0] = w0;
  ow[1] = w1;
  ow[2] = w2;
  DirectANNBasisSet dabs(wghts);
  DirectANNModel dam(wghts,ow);
  VecDbl x(2);
  x[0] = x0;
  x[1] = x1;
  cout << "eval: " << dam(x) << endl;
  VecDbl grad = dam.gradient(x);
  cout << "ddam/d0: " << grad[0] << endl;
}

const unsigned GRIDPOINTS = 50;
void SurfpackModelTest::generalDerivativeTest(const SurfpackModel& model, const AxesBounds& ab)
{
  // Randomly sample a few  points
  // The derivative of the model will be evaluated at each of these
  // points.  A tangent plane will be constructed at each point, based
  // on the gradient.  Plots will show projections of the model and the
  // tangent plane along each dimension
  const unsigned nsamples = 1;
  vector< string > test_functions;
  SurfData* randsd = 0;
  SurfpackInterface::CreateSample(randsd,ab,nsamples,test_functions);
  

  // Prepare a basis set for a hyperplane model 
  // This will be the tangent (hyper)plane: 
  //   y = c0 + c1*x1 + c2*x2 + ... + cn*xn;
  LRMBasisSet hpbs; // hyperplane basis set
  VecUns list(1);
  for (unsigned i = 0; i < randsd->xSize(); i++) {
    list[0] = i;
    hpbs.bases.push_back(list);
  }
  // Put the unity basis on the end so that the coefficent for that
  // basis can just be tacked on the back of the gradient
  hpbs.bases.push_back(VecUns()); // Unity

  // For each of the randomly sampled points
  for (unsigned i = 0; i < randsd->size(); i++) {
    VecDbl pt = (*randsd)(i);
    // Compute the model value and gradient at that point
    double resp = model(pt);
    VecDbl grad = model.gradient(pt);
    // solve for b in y = m1x1 + m2x2 + ... + b;
    // The gradient of the function gives the slope of the tangent
    // plane.  We just need to find the appropriate offset from the origin.
    double b = resp;
    for (unsigned j = 0; j < pt.size(); j++) {
      b -= pt[j]*grad[j];
    }

    // Now create a tangent hyperplane to the surface
    VecDbl hyp_coeff = grad;
    hyp_coeff.push_back(b);
    LinearRegressionModel tang_plane(pt.size(),hpbs,hyp_coeff);

    // Now prepare a set of axes, so that we can plot a slice of both
    // functions along each dimension
    // This is a set of axes along which each dimension is fixed
    vector<AxesBounds::Axis> axes(ab.size());
    for (unsigned j = 0; j < pt.size(); j++) {
      axes[j].min = axes[j].max = pt[j];
    }
    
    for (unsigned dim = 0; dim < pt.size(); dim++) {
      vector<AxesBounds::Axis> var_axes = axes;
      var_axes[dim].min = ab[dim].min;
      var_axes[dim].max = ab[dim].max;
      AxesBounds var_ab(var_axes);
      // Number of points per dim will be 1 for all but the active dimension
      VecUns var_dims(ab.size(),1);
      var_dims[dim] = GRIDPOINTS;

      // Now create the data set that varies only along the active dimension
      SurfData* test_data = 0;
      vector<string> test_functions;
      SurfpackInterface::CreateSample(test_data,var_ab,var_dims,test_functions);
      SurfData& td = *test_data;
      // Now evaluate both models on this data
      ///\\todo Creat a base model class that can do operator()(SurfData&)
      VecDbl model_resp(td.size());
      VecDbl hyp_resp(td.size());
      for (unsigned k = 0; k < td.size(); k++) {
        model_resp[k] = model(td(k));
        hyp_resp[k] = tang_plane(td(k));
      }

      // Now generate all the output files
      unsigned y1 = td.addResponse(model_resp) + td.xSize() + 1; // gnuplot counts cols from 1
      unsigned y2 = td.addResponse(hyp_resp) + td.xSize() + 1;
      unsigned xdim = dim + 1;
      string s1("deriv");
      ostringstream os;
      os << s1 << i << "_" << dim ;
      string plotfilename = os.str();
      // This file will contain the test data (which varies along a single dim),
      // model response, and values of the tangent plane
      string datafilename = plotfilename + ".spd";
      // This file contains only the randomly generated point at which
      // the gradient was evaluated (to be used for emphasis in plot)
      string ptfilename = plotfilename + "pt";
      std::ofstream ptfile(ptfilename.c_str(),ios::out);
      ptfile << pt[dim] << " " << resp << endl;
      ptfile.close(); 
      td.write(datafilename);
      generateThePlot(plotfilename,plotfilename,xdim,y1,y2); 
      delete test_data;
    }
    
  } 
  delete randsd;
}

void SurfpackModelTest::generateThePlot(const std::string& plotname, const std::string& datafilename, int x, int y1, int y2)
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

