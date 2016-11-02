#include "surfpack_system_headers.h"
#include "surfpack.h"
#include "LinearRegressionModel.h"
#include "SurfData.h"
#include "ModelScaler.h"

using std::cout;
using std::endl;
using std::accumulate;
using std::string;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(LinearRegressionModel)
#endif


double LRMBasisSet::eval(unsigned index, const VecDbl& x) const
{
  assert(index < bases.size());
  double result = 1.0;
  for(VecUnsIt it = bases[index].begin(); it != bases[index].end(); ++it) {
    if (*it >= x.size()) cout << index << " " << *it << endl;
    assert(*it < x.size());
    result *= x[*it];
  } 
  return result;
}

double LRMBasisSet::deriv(unsigned index, const VecDbl& x, const VecUns& vars) const
{
  VecUns counts(x.size(),0);
  for(VecUnsIt it = bases[index].begin(); it != bases[index].end(); ++it) {
    assert(*it < x.size());
    counts[*it]++;
  } 
  double coeff = 1.0;
  for(VecUns::const_iterator it = vars.begin(); it != vars.end(); ++it) {
    assert(*it < x.size());
    unsigned index = *it;
    // Taken derivative with respect to this variable too many times
    if (!counts[index]) return 0.0; 
    coeff *= counts[index]--;
  }
  unsigned sum = accumulate(counts.begin(),counts.end(),0);
  if (sum == 0) return coeff; // no vars left
  double term = 1.0;
  for(unsigned v = 0; v < counts.size(); v++) {
    for (unsigned c = 0; c < counts[v]; c++) {
      term *= x[v];
    }
  }
  //cout << "Deriv: " << index << " X: ";
  //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout," "));
  //cout << " Vars: " ;
  //copy(vars.begin(),vars.end(),std::ostream_iterator<double>(cout," "));
  //cout << " Val: " << coeff*term << endl;
  return coeff*term;
}
  
std::string LRMBasisSet::asString() const
{
  std::ostringstream os;
  for(VecVecUns::const_iterator it = bases.begin(); it != bases.end(); ++it) {
    if (it->empty()) { os << "Unity\n"; continue; }
    copy(it->begin(),it->end(),std::ostream_iterator<unsigned>(os," "));
    os << "\n";
  }
  return os.str();
}

void LRMBasisSet::add(const std::string& s_basis)
{
  bases.push_back(surfpack::toVec<unsigned>(s_basis));
}


LinearRegressionModel::LinearRegressionModel(const unsigned dims, 
  const LRMBasisSet& bs_in, const VecDbl& coeffs_in)
  : SurfpackModel(dims), bs(bs_in), coeffs(coeffs_in)
{
  assert(bs.bases.size() == coeffs.size());
}

double LinearRegressionModel::evaluate(const VecDbl& x) const
{
  static int times_called = 0;
  assert(coeffs.size() == bs.bases.size());
  double sum = 0;
  for (unsigned i = 0; i < coeffs.size(); i++) {
    sum += coeffs[i]*bs.eval(i,x);
  }
  //cout << "LinearRegression: times called: " << ++times_called << endl;
  return sum;
}

VecDbl LinearRegressionModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  cout << "IN gradient x[0] = " << x[0] << endl;
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

std::string LinearRegressionModel::asString() const
{
  std::ostringstream os;
  unsigned num_bases = bs.size();
  unsigned num_vars = ndims;
  os << "-----\n";
  os << "Surfpack polynomial model\n";
  os << "f(x) = sum_k{c_k * prod_k[x(i) ^ p(k,i)]}; where\n";
  os << "\ninputs = " << num_vars << "\n";
  os << "bases = " << num_bases << "\n";
  os << "\nc (1 x bases) =\n";
  os << std::scientific << std::setprecision(16);
  for(unsigned i=0; i<num_bases; ++i)
    os << std::setw(23) << coeffs[i] << " ";
  os << "\n\np (bases x inputs) = \n";
  os << std::fixed << std::setprecision(0);
  for(VecVecUns::const_iterator it = bs.bases.begin(); it != bs.bases.end(); ++it) {
    for(unsigned i = 0; i < num_vars; i++)
       os << std::setw(3) << std::count(it->begin(), it->end(), i) << " ";
    os << "\n";
  }
  os << "-----\n";
  return os.str();
}


///////////////////////////////////////////////////////////
///	Polynomial (LinearRegression) Model Factory
///////////////////////////////////////////////////////////

VecDbl LinearRegressionModelFactory::lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd)
{
  MtxDbl A(ssd.size(),bs.size(),true);
  for (unsigned i = 0; i < ssd.size(); i++) {
    for (unsigned j = 0; j < bs.size(); j++) {
      A(i,j) = bs.eval(j,ssd(i));
    }
  }
  VecDbl b = ssd.getResponses();
  VecDbl x(bs.size());
  if (eqConRHS.empty()) {
    surfpack::linearSystemLeastSquares(A,x,b);
  } else {
    surfpack::leastSquaresWithEqualityConstraints(A,x,b,eqConLHS,eqConRHS);
  }
  return x; 
}

LRMBasisSet LinearRegressionModelFactory::CreateLRM(unsigned order, 
  unsigned dims)
{
  LRMBasisSet bs;
  bs.add(std::string(""));
  std::deque<Term> q;
  q.push_front(Term(VecUns()));
  while (!q.empty()) {
    Term& t = q.front();
    VecUns& v = t.vars;
    if (v.size() < order && !t.color) { // extendable
      t.color = true; // only extend it once
      Term new_term = Term(v);
      if (v.empty()) new_term.vars.push_back(0);
      else new_term.vars.push_back(v.back());
      bs.bases.push_back(new_term.vars);
      q.push_front(new_term);
    } else if (!v.empty() && v.back() < dims-1) {
      v.back()++;
      bs.bases.push_back(v);
      t.color = false;
    } else {
      q.pop_front();
    }
  }
  return bs;
} 


SurfpackModel* LinearRegressionModelFactory::Create(const SurfData& sd)
{
  // historically the constraintPoint was not scaled and equality
  // constraints were set before Build (therefore Create) was called,
  // so we call here before other operations
  // TODO: sort out scaling and make consistent across models
  setEqualityConstraints(sd.getConstraintPoint());

  //ModelScaler* ms = NormalizingScaler::Create(sd);
  ModelScaler* ms = NonScaler::Create(sd);
  ScaledSurfData ssd(*ms,sd);
  
  LRMBasisSet bs = CreateLRM(order,sd.xSize());
  //cout << bs.asString() << endl;
  //cout << "sd size: " << sd.size() << " bs size: " << bs.size() <<  endl;
  VecDbl coeffs = lrmSolve(bs,ssd);
  //copy(coeffs.begin(),coeffs.end(),std::ostream_iterator<double>(cout,"|"));
  //cout << "\n";
  SurfpackModel* lrm = new LinearRegressionModel(sd.xSize(),bs,coeffs);
  lrm->scaler(ms);
  delete ms;
  return lrm;
}

unsigned LinearRegressionModelFactory::minPointsRequired()
{
  config();
  LRMBasisSet bs = CreateLRM(order,ndims);
  return bs.size();
}

unsigned LinearRegressionModelFactory::recommendedNumPoints()
{
  return LinearRegressionModelFactory::minPointsRequired();
}

bool LinearRegressionModelFactory::supports_constraints()
{
  return true;
}

/// total data from points and constraint must be sufficient to build
void LinearRegressionModelFactory::sufficient_data(const SurfData& sd)
{
  if (sd.size() + sd.numConstraints() < minPointsRequired()) {
    std::ostringstream not_enough;
    not_enough << "Not enough Points: ";
    not_enough << "size of data = " << sd.size();
    not_enough << ", size of constraints data = " << sd.numConstraints();
    not_enough << ", minPointsRequired = " << minPointsRequired();
    throw(not_enough.str());
  }
}


LinearRegressionModelFactory::LinearRegressionModelFactory()
  : SurfpackModelFactory(), order(2)
{

}

LinearRegressionModelFactory::LinearRegressionModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), order(2)
{

}

void LinearRegressionModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["order"];
  if (strarg != "") order = std::atoi(strarg.c_str());
}

void LinearRegressionModelFactory::setEqualityConstraints(const SurfPoint& sp)
{
  unsigned short asv = 0;
  if (sp.fSize() > 0) asv |= 1;
  if (sp.fGradientsSize() > 0) asv |= 2;
  if (sp.fHessiansSize() > 0) asv |= 4;
  if (asv == 0) return; // There are no constraints

  config();
  LRMBasisSet bs = CreateLRM(order,ndims);
  VecDbl coefficients(bs.size());

  unsigned numConstraints = 0;
  if (asv & 1) numConstraints += 1; // value at a particular point
  if (asv & 2) numConstraints += ndims; // gradient at a point
  if (asv & 4) numConstraints += (ndims*ndims+ndims)/2; // hessian at a point
  eqConRHS.resize( numConstraints );
  // Must compute number of terms first
  //MtxDbl temp(eqConRHS.size(),coefficients.size(),true);
  //eqConLHS = temp;
  eqConLHS.reshape(eqConRHS.size(),coefficients.size());
  // Marks the index of the next constraint to be added (necessary since
  // indices of e.g. the gradient constraints will be different depending on
  // whether or not the value constraint is used
  unsigned index = 0;
  // If requested, add the equality constraint for the point value
  if (asv & 1) {
    for (unsigned i = 0; i < bs.size(); i++) {
      eqConLHS(index,i) = bs.eval(i,sp.X());
    }
    eqConRHS[index] = sp.F();
    ++index;
  }

  // If requested, add the equality constraints for the gradient
  if (asv & 2) {
    const VecDbl& gradient = sp.fGradient();
    assert(gradient.size() == ndims);
    VecUns factorCounts;
    VecUns diff_vars(1); // Holds index of var to differentiate w.r.t.
    for (unsigned dif_var = 0; dif_var < ndims; dif_var++ ) {
      diff_vars[0] = dif_var;
      for (unsigned i = 0; i < bs.size(); i++) {
        eqConLHS(index,i) = bs.deriv(i,sp.X(), diff_vars);
      }
      eqConRHS[index] = gradient[dif_var];
      ++index;
    }
  }

  // If requested, add the equality constraints for the hessian
  if (asv & 4) {
    const MtxDbl& hessian = sp.fHessian();
    assert(hessian.getNCols() == ndims);
    assert(hessian.getNRows() == ndims);
    VecUns factorCounts;
    VecUns diff_vars(2); // Holds indices of vars to differentiate w.r.t.
    for (unsigned dif_var1 = 0; dif_var1 < ndims; dif_var1++ ) {
      diff_vars[0] = dif_var1;
      for (unsigned dif_var2 = dif_var1; dif_var2 < ndims; dif_var2++ ) {
        diff_vars[1] = dif_var2;
        for (unsigned i = 0; i < bs.size(); i++) {
          eqConLHS(index,i) = bs.deriv(i,sp.X(), diff_vars);
        }
        eqConRHS[index] = hessian(dif_var1,dif_var2);
        ++index;
      } // dif_var2
    } // dif_var1
  } // if hessian needed
}
