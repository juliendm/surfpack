#include "surfpack.h"
#include "MovingLeastSquaresModel.h"

using std::cout;
using std::endl;
using std::string;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(MovingLeastSquaresModel)
#endif


double weight(const VecDbl xi, const VecDbl x, unsigned continuity = 1, 
	      double radius = 1.0)
{
  assert(continuity > 0);
  assert(continuity < 4);
  
  double weight = 0.0;
  double rho = surfpack::euclideanDistance(xi,x)/radius;

  switch (continuity) {

  case 1:
    weight = exp(-1*(pow(rho,2.0)))/(pow(rho,2.0)+0.001);
    // weight = (rho > 1.0) ? 0.0 : 1.0-3.0*rho*rho+2.0*pow(rho,3.0);
    break;

  case 2:
    weight = (rho > 1.0) ? 0.0
      : 1.0-10.0*pow(rho,3.0)+15.0*pow(rho,4.0)-6.0*pow(rho,5.0);
    break;

  case 3:
    weight = (rho > 1.0) ? 0.0 
      : 1.0-35.0*pow(rho,4.0)+84.0*pow(rho,5.0)-70.0*pow(rho,6.0)+20.0*pow(rho,7.0);
    break;

  }

  return weight;
}

MovingLeastSquaresModel::MovingLeastSquaresModel(const SurfData& sd_in, const LRMBasisSet& bs_in, unsigned continuity_in)
  : SurfpackModel(sd_in.xSize()), sd(sd_in), bs(bs_in), continuity(continuity_in)
{
  assert(continuity > 0);
  assert(continuity < 4);
}

double MovingLeastSquaresModel::evaluate(const VecDbl& x) const
{
  unsigned nbases = bs.bases.size();
  MtxDbl A(nbases,nbases,true);
  VecDbl By(nbases,0.0); // B(x) = Pt(x)*w(x); By = B(x)*y;
  VecDbl resps = sd.getResponses();
  // # of data points must be at least as great as the number of basis functions
  assert(resps.size() >= bs.size()); 
  for (unsigned i = 0; i < nbases; i++) {
    for (unsigned j = 0; j < nbases; j++) {
      A(i,j) = 0.0;
      for (unsigned k = 0; k < sd.size(); k++) {
        A(i,j) += bs.eval(i,sd(k))*bs.eval(j,sd(k))*weight(sd(k),x,continuity);
        if (!j) By[i] += bs.eval(i,sd(k))*weight(sd(k),x,continuity)*resps[k];
      }
    }
  }
  //  cout << "Amatrix" << "\n";  
  //for (unsigned i = 0; i < nbases; i++) {
  //  for (unsigned j = 0; j < nbases; j++) {
  //    cout << A(i,j) << " " ;
  //  }
  //  cout << "\n";
  //} 

  //cout << "Bmatrix" << "\n";  
  //for (unsigned i = 0; i < nbases; i++) 
  // cout << By[i] << " " ;
  //cout << "\n";

  //VecDbl temp;
  surfpack::linearSystemLeastSquares(A,coeffs,By);
  //coeffs = temp;
  double sum = 0.0;
  for(unsigned i = 0; i < nbases; i++) {
    sum += bs.eval(i,x)*coeffs[i];
    //cout << "Coefficients " << i << "= " << coeffs[i];
    //cout << "basis " << i << "= " << bs.eval(i,x);
  }
  return sum;
}

/// Currently set up so that operator() must be called immediately before
/// Not good assumption
///\todo Be able to recompute coefficients so that prior op() call is not req.
VecDbl MovingLeastSquaresModel::gradient(const VecDbl& x) const
{
  ///To get the coefficients for the derivatives of the basis functions
  ///we need to call the evaluate function.  This is potentially inefficient,
  ///since it is likely that the evaluate function would have been called
  ///prior to this call.  But it should at least be correct.
  ///\todo Remove inefficiency of having to call evaluate when gradient is called.
  (*this)(x);
  /// code copied straight from LRM
  assert(!x.empty());
  assert(coeffs.size() == bs.bases.size());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl result(x.size(),0.0);
  for (unsigned i = 0; i < x.size(); i++) {
    diff_var[0] = i;
    for (unsigned j = 0; j < bs.bases.size(); j++) {
      result[i] += coeffs[j]*bs.deriv(j,x,diff_var);
    }
  }
  return result;
}

std::string MovingLeastSquaresModel::asString() const
{
  std::ostringstream os;
  os << "\nbases:\n" << bs.asString() << "\n";
  os << "\ncontinuity: " << this->continuity << endl;
  return os.str();
}

///////////////////////////////////////////////////////////
///	Moving Least Squares Model Factory
///////////////////////////////////////////////////////////

MovingLeastSquaresModelFactory::MovingLeastSquaresModelFactory()
  : SurfpackModelFactory(), weight(1), order(2)
{

}

MovingLeastSquaresModelFactory::MovingLeastSquaresModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), weight(1), order(2)
{

}

void MovingLeastSquaresModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["weight"];
  if (strarg != "") weight = std::atoi(strarg.c_str());
  strarg = params["order"];
  if (strarg != "") order = std::atoi(strarg.c_str());
}

SurfpackModel* MovingLeastSquaresModelFactory::Create(const SurfData& sd)
{
  LRMBasisSet lrmbs = LinearRegressionModelFactory::CreateLRM(order,ndims);
  SurfpackModel* sm = new MovingLeastSquaresModel(sd, lrmbs, weight); 
  assert(sm);
  return sm; 
}
