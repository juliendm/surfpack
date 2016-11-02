#include "surfpack_system_headers.h"
#include "surfpack.h"
#include "MarsModel.h"
#include "SurfData.h"

#ifdef HAVE_CONFIG_H
#define FMODM_F77 F77_FUNC(fmodm,FMODM)
#define MARS_F77  F77_FUNC(mars,MARS)
#else
#include "surf77_config.h"
#define FMODM_F77 SURF77_GLOBAL(fmodm,FMODM)
#define MARS_F77  SURF77_GLOBAL(mars,MARS)
#endif

extern "C" { // prevent C++ name mangling
  void FMODM_F77(int&, int&, real*, real*, int*, real*, real*);

  void MARS_F77(int&, int&, real*, real*, real*, int&, int&, int*,
                real*, int*, real*, double*, int*);
} // end extern "C"


using std::cout;
using std::endl;
using std::string;
using std::max;
using std::memcpy;
using std::memset;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(MarsModel)
#endif


MarsModel::MarsModel(const unsigned dims, real* fm_in, int fmsize, int* im_in, 
  int imsize, int interp)
  : SurfpackModel(dims), fm(fmsize), im(imsize), interpolation(interp)
{
  memcpy(&fm[0],fm_in,sizeof(real)*fmsize);
  memcpy(&im[0],im_in,sizeof(int)*imsize);
}

double MarsModel::evaluate(const VecDbl& x) const
{
  int nval = 1;
  real* xVector = new real[static_cast<int>(x.size())];
  for (int i = 0; i < static_cast<int>(x.size()); i++) {
    xVector[i] = static_cast<real>(x[i]);
  }
  real* sp = new real[2];
  real* f = new real[1];
  // (BMA, 4/17/2007): The following explicit initializations added due to
  // similar need with fm and im in build(...) below.  Conservative.
  memset(sp, 0, 2*sizeof(real));
  f[0] = 0;
  int continuity_level = interpolation;
  FMODM_F77(continuity_level,nval,xVector,const_cast<real*>(&fm[0]),
    const_cast<int*>(&im[0]),f,sp);
  delete [] sp;
  delete [] xVector;
  real result = *f;
  delete [] f;
  return result;
}

VecDbl MarsModel::gradient(const VecDbl& x) const
{
  throw string("Mars does not currently provide analytical gradients");
}

std::string MarsModel::asString() const
{
  std::ostringstream os;
  os << "Mars model\n" ;
  return os.str();
}



///////////////////////////////////////////////////////////
///	MARS Model Factory
///////////////////////////////////////////////////////////

SurfpackModel* MarsModelFactory::Create(const SurfData& sd)
{
  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();
  delete [] xMatrix;
  delete [] fm;
  delete [] im;
  int nmcv = 0;
  int ntcv = 0;
  //int max_bases = 15;
  //int max_interactions = 2;
  n = static_cast<int>(sd.size());
  np = static_cast<int>(sd.xSize());
  xMatrix = new real[n*np];
  real* y = new real[n];
  real* w = new real[n];
  int* lx = new int[np];
  int fmsize = 3+max_bases*(5*max_interactions+nmcv+6)+2*np+ntcv;
  int imsize = 21+max_bases*(3*max_interactions+8);
  fm = new real[fmsize];
  im = new int[imsize];
  real* sp = new real[2*(n*(max(max_bases+1,2)+3)+
			 max(3*n+5*max_bases+np,max(2*np,4*n)))+
		      2*np+4*max_bases];
  double* dp = new double[2*(max(n*max_bases,(max_bases+1)*(max_bases+1))+
			     max((max_bases+2)*(nmcv+3),4*max_bases))];
  int* mm = new int[2*(n*np+2*max(max_interactions,nmcv))];

  // (BMA, 04/17/2007): The following explicit initializations added due to
  // valgrind reports of accessing uninitialized memory in fm.  Perhaps due to
  // old FORTRAN compiler convention of zeroing variables.  This is likely more
  // conservative than needed, but should be safer (at a cost).
  memset( fm, 0, (3+max_bases*(5*max_interactions+nmcv+6)+2*np+ntcv)*
	         sizeof(real) );
  memset( im, 0, (21+max_bases*(3*max_interactions+8))*sizeof(int) );
  memset( sp, 0, (2*(n*(max(max_bases+1,2)+3) + 
		     max(3*n+5*max_bases+np,max(2*np,4*n))) +
		  2*np+4*max_bases) * sizeof(real) );
  memset( dp, 0, 2*(max(n*max_bases,(max_bases+1)*(max_bases+1))+
		    max((max_bases+2)*(nmcv+3),4*max_bases))*sizeof(double) );
  memset( mm, 0, 2*(n*np+2*max(max_interactions,nmcv))*sizeof(int) );

  //unsigned pts = data.size();
  //for (unsigned i = 0; i < pts; i++) {
  for (int i = 0; i < n; i++) {
    //cout << current->getF(responseIndex) << endl;
    //const vector<double> domain = current.X();
    for (int j = 0; j < np; j++) {
      xMatrix[j*n+i] = static_cast<real>(sd(i,j)); 
    }
    y[i] = static_cast<real>(sd.getResponse(i));
    w[i] = 1.0f;
  } 
  // Specify each variable to be 'unrestricted'
  for (int k = 0; k < np; k++) {
    lx[k] = 1;
  }
  //printMatrix(xMatrix,n,np,cout);
  //printMatrix(w,n,1,cout);
  //printMatrix(y,n,1,cout);
  //printIntMatrix(lx,np,1,cout);
  MARS_F77(n,np,xMatrix,y,w,max_bases,max_interactions,lx,fm,im,sp,dp,mm);
  SurfpackModel* model = new MarsModel(ndims,fm,fmsize,im,imsize,interpolation);
  delete [] y;
  delete [] w;
  delete [] lx;
  delete [] sp;
  delete [] dp;
  delete [] mm;
  assert(model);
  return model;
}

MarsModelFactory::MarsModelFactory()
  : SurfpackModelFactory(),
  xMatrix(0), fm(0), im(0), max_interactions(2), max_bases(25)
{

}

MarsModelFactory::MarsModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args),
  xMatrix(0), fm(0), im(0), max_interactions(2), max_bases(25)
{

}

void MarsModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["max_bases"];
  if (strarg != "") max_bases = std::atoi(strarg.c_str());
  strarg = params["max_interactions"];
  if (strarg != "") max_interactions = std::atoi(strarg.c_str());
  strarg = params["interpolation"];
  if (strarg == "linear") interpolation = 1;
  else if (strarg == "cubic") interpolation = 2;
  else if (strarg != "") throw string("Mars interpolation must be linear or cubic");
}
