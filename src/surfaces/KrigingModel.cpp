#include "surfpack_system_headers.h"
#include "KrigingModel.h"
#include "SurfpackMatrix.h"
#include "SurfData.h"
#include "surfpack.h"
#include "ModelScaler.h"
#include "AxesBounds.h"
#include "ModelFitness.h"

using std::cout;
using std::endl;
using std::copy;
using std::vector;
using std::string;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(KrigingModel)
#endif


void KrigingModel::
surfdata_to_nkm_surfdata(const SurfData& sd, nkm::SurfData& nkm_sd)
{
  unsigned num_points = sd.size();
  unsigned x_size = sd.xSize();
  // BMA: NKM allows the data to have multiple responses, but the
  // Surfpack architecture (and NKM) assumes the model is built over a
  // single active response index.  Here, we only map in the data for
  // a single response instead of using the nkm::SurfData::setIOut()
  // feature, to avoid duplicating the response data for every model.
  // We leave the general loops provided by KRD until future
  // improvements.
  unsigned f_active = sd.getDefaultIndex();
  unsigned f_size = 1;
  //  unsigned f_size = sd.fSize();

  nkm::MtxDbl XR(x_size,num_points), Y(f_size, num_points);
  std::vector<std::vector<nkm::MtxDbl> > derY(f_size);
  // track the active derivative order for each function
  nkm::MtxInt der_order(f_size,1); 
  der_order.zero(); //set contents to zero

  if(num_points>0) {
    //could increment der_order independently for each dimension but old surfdata/surfpoint does not support this, nkm::SurfData supports arbitrarily high order derivatives.

    const SurfPoint& point=sd[0];
    if(point.fGradientsSize() > 0) {
      for(unsigned f_index=0; f_index<f_size; ++f_index)
	++der_order(f_index,0);
      if(point.fHessiansSize() > 0)
	for(unsigned f_index=0; f_index<f_size; ++f_index)
	  ++der_order(f_index,0);
    }
    
    for(unsigned f_index=0; f_index<f_size; ++f_index) {
      derY[f_index].resize(der_order(f_index,0)+1);
      for(unsigned der_order_index=1; static_cast<int>(der_order_index)<=der_order(f_index,0); ++der_order_index)
	derY[f_index][der_order_index].newSize(nkm::num_multi_dim_poly_coef(x_size,-der_order_index),num_points);
    }
  }


  for (unsigned point_index=0; point_index<num_points; ++point_index) {
    const SurfPoint& point = sd[point_index];
    const VecDbl x = point.X();
    for (unsigned x_index=0; x_index<x_size; ++x_index)
      XR(x_index,point_index) = x[x_index];
    for (unsigned f_index=0; f_index<f_size; ++f_index) 
      Y(f_index,point_index) = point.F(f_active);


    // given NKM ordering of derY, probably need to loop differently
    // to populate the matrices; this just for demo

    // example of mapping first derivatives
    // there should be 0 or f_size gradients (could throw error)
    if (point.fGradientsSize() > 0) {
      for (unsigned f_index=0; f_index < f_size; ++f_index) {
	assert(der_order(f_index,0)>=1);  //could change this to a throw
	const vector<double>& sd_gradient = point.fGradient(f_active);
	//cout << "Surfpack gradient for point " << point_index << ", function "
	//     << f_index << ": [ ";
	for (unsigned x_index=0; x_index < x_size; ++x_index) {
	  // accessing each gradient element
	  //cout << sd_gradient[x_index] << " ";
	  derY[f_index][1](x_index,point_index) =sd_gradient[x_index];
	}
	//cout << "]" << endl;
      }
    }
    else{
      for (unsigned f_index=0; f_index < f_size; ++f_index) 
	assert(der_order(f_index,0)==0);  //could change this to a throw
    }

    // example of mapping second derivatives
    // there should be 0 or f_size Hessians (could throw error)
    if (point.fHessiansSize() > 0) {
      for (unsigned f_index=0; f_index<f_size; ++f_index) {
	assert(der_order(f_index,0)>=2);  //could change this to a throw

	const SurfpackMatrix<double>& sd_hessian = point.fHessian(f_active);
	//cout << "Surfpack Hessian for point " << point_index << ", function "
	//     << f_index << " is:\n";
	unsigned der_index=0;
	for (unsigned xj_index=0; xj_index < x_size; ++xj_index) 
	  for (unsigned xk_index=xj_index; xk_index < x_size; ++xk_index, ++der_index)
	    derY[f_index][2](der_index,point_index)=
	      sd_hessian(xj_index,xk_index);

	//for (unsigned xj_index=0; xj_index < x_size; ++xj_index) {
	//  for (unsigned xk_index=0; xk_index < x_size; ++xk_index) {
	//    // accessing each Hessian element
	//    cout << sd_hessian(xj_index, xk_index) << " ";
	//  }
	//  cout << "\n";
	//}
	//cout << endl;
      }
    }
    else{
      for(unsigned f_index=0; f_index<f_size; ++f_index)
	assert(der_order(f_index,0)<2);
    }
  }

  // If present, add constraintPoint as well, indicating index in
  // arguments as needed, adding to f, grad, Hess
  const SurfPoint& constraint_point = sd.getConstraintPoint();

  nkm_sd = nkm::SurfData(XR, Y, der_order, derY);

  // Could use this, but instead opting to not copy needless columns
  // of data into each model...
  //  nkm_sd.setIOut(f_active_index);
}


KrigingModel::KrigingModel(const SurfData& sd, const ParamMap& args)
  : SurfpackModel(sd.xSize()), nkmKrigingModel(NULL)
{
  nkm::SurfData nkmSurfData;
  surfdata_to_nkm_surfdata(sd, nkmSurfData);
  /*
  int der_order=0;
  ParamMap::const_iterator param_it;
  param_it = args.find("derivative_order");
  if (param_it != args.end() && param_it->second.size() > 0) {
    der_order=std::atoi(param_it->second.c_str());
  }
  if(der_order<0) {
    std::cerr << "error in KrigingModel::KrigingModel()\nyou specified a negative derivative order, but derivative_order must be a non-negative integer" << std::endl;
    assert(0<=der_order);
  }
  if(der_order>1) {
    std::cerr << "error in KrigingModel::KrigingModel()\nyou specified that you want to build the Kriging Model with derivative order greater than 1.  Currently only regular kriging and gradient enhacned Kriging is implemented" << std::endl;
    assert(der_order<=1);
  }
  if(der_order==0)
  else if(der_order==1) {
    if(nkmSurfData.getDerOrder()<1) {
      std::cerr << "error in KrigingModel::KrigingModel()\nyou requested that we build a Gradient Enhanced Kriging Model but did not provide gradient information" << std::endl;
      assert(1<=nkmSurfData.getDerOrder());
    }
    nkmKrigingModel = new nkm::GradKrigingModel(nkmSurfData, args);
    //std::cerr << "built Gradient Enhanced Kriging Model" << std::endl;
  }
  */
  nkmKrigingModel = new nkm::KrigingModel(nkmSurfData, args);
  nkmKrigingModel->create();

}


KrigingModel::~KrigingModel()
{
  if (nkmKrigingModel)
    delete nkmKrigingModel;
}


double KrigingModel::evaluate(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(ndims,1);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(i,0) = x[i];
  return (nkmKrigingModel->evaluate(nkm_x));
}


double KrigingModel::variance(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(ndims,1);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(i,0) = x[i];
  //double adj_var=nkmKrigingModel->eval_variance(nkm_x);
  //printf("KrigingModel::variance() adj_var=%g\n",adj_var);
  //return (adj_var);
  return (nkmKrigingModel->eval_variance(nkm_x));
}


VecDbl KrigingModel::gradient(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(ndims,1);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(i,0) = x[i];

  nkm::MtxDbl nkm_d1y(ndims,0);
  nkmKrigingModel->evaluate_d1y(nkm_d1y, nkm_x);

  VecDbl d1y(ndims, 0.0); 
  for(size_t i=0; i<ndims; ++i)
    d1y[i] = nkm_d1y(i,0);

  return d1y;
}

MtxDbl KrigingModel::hessian(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(ndims,1);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(i,0) = x[i];

  int num_lower_elem=(ndims+1)*ndims/2;
  nkm::MtxDbl nkm_d2y(num_lower_elem,1);
  nkmKrigingModel->evaluate_d2y(nkm_d2y, nkm_x);

  MtxDbl d2y(ndims, ndims, 0.0); 
  int k=0;
  for(size_t j=0; j<ndims; ++j) {
    d2y(j,j)=nkm_d2y(k,0);
    k++;
    for(size_t i=j+1; i<ndims; ++i) {
      d2y(i,j) = nkm_d2y(k,0);
      d2y(j,i) = d2y(i,j);
      k++;
    }
  }

  return d2y;
}


std::string KrigingModel::asString() const
{

  return (nkmKrigingModel->asString());
  // TODO: be able to write NKM_KrigingModel as a string
  //assert(false);


}


//KrigingModel KrigingModel::Create(const SurfData& sd)
//{
//  ConminKriging ck(sd);
//  VecDbl best_guess(sd.xSize(),1.0);
//  unsigned max_iter = 100; 
//  double opt_val;
//  VecDbl rhs;
//  ck.optimize(best_guess,opt_val,max_iter); 
//  KrigingBasisSet kbs(SurfData::asVecVecDbl(sd),best_guess);
//}

///////////////////////////////////////////////////////////
/// 	Kriging Model Factory	
///////////////////////////////////////////////////////////

KrigingModelFactory::KrigingModelFactory()
  : SurfpackModelFactory()
{

}

KrigingModelFactory::KrigingModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args)
{

}

void KrigingModelFactory::config()
{
  SurfpackModelFactory::config();
}

bool KrigingModelFactory::supports_constraints()
{
  // TODO: Kriging model will soon map the anchor point
  //  return true;
  return false;
}

void KrigingModelFactory::sufficient_data(const SurfData& sd)
{
  // NKM manages error checking, so this always passes
  return;
}


typedef std::pair<double,VecDbl> KMPair;
SurfpackModel* KrigingModelFactory::Create(const SurfData& sd)
{


  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();

  return new KrigingModel(sd, params);

}
