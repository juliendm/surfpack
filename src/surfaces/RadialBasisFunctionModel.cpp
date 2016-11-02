#include "surfpack_system_headers.h"
#include "RadialBasisFunctionModel.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "ModelFitness.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::max;
using std::min;
using surfpack::shared_rng;
using surfpack::fromVec;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(RadialBasisFunctionModel)
#endif


SurfPoint computeCentroid(const SurfData& sd)
{
  assert(sd.size());
  assert(sd.xSize());
  VecDbl center(sd.xSize(),0.0);
  for (unsigned pt = 0; pt < sd.size(); pt++) {
    for (unsigned dim = 0; dim < sd.xSize(); dim++) {
      center[dim] += sd(pt,dim);
    }
  }
  for (unsigned dim = 0; dim < center.size(); dim++) {
    center[dim] /= sd.size();
  }
  return SurfPoint(center);
}

void updateCentroid(VecDbl& centroid, const VecDbl& newpt, unsigned weight)
{
  assert(centroid.size() == newpt.size());
  for (unsigned i = 0; i < centroid.size(); i++) {
    if (!weight) centroid[i] = newpt[i];
    else centroid[i] = (weight*centroid[i]+newpt[i])/(weight+1);
  }
}

unsigned findClosest(const SurfData& sd, VecDbl pt)
{
  assert(sd.size());
  double mindist = surfpack::euclideanDistance(sd(0),pt);
  unsigned argmin = 0;
  for (unsigned i = 1; i < sd.size(); i++) {
    double distance = surfpack::euclideanDistance(sd(i),pt);
    if (distance < mindist) {
      mindist = distance;
      argmin = i;
    }
  }
  return argmin;
}

SurfData radii(const SurfData& generators)
{
  SurfData result;
  for (unsigned i = 0; i < generators.size(); i++) {
    VecDbl radius(generators.xSize(),std::numeric_limits<double>::max());
    for (unsigned j = 0; j < generators.size(); j++) {
      if (i != j) {
        for (unsigned dim = 0; dim < generators.xSize(); dim++) {
          double distance = fabs(generators(i,dim)-generators(j,dim));
          if (distance < radius[dim]) radius[dim] = distance;
        }
      }
    }
    result.addPoint(SurfPoint(radius));
  }
  return result;
}

SurfData cvts(const AxesBounds& ab, unsigned ngenerators, unsigned ninfluencers,
  double minalpha = .5, double maxalpha = .99)
{
  assert(ninfluencers > ngenerators);
  SurfData* generators = ab.sampleMonteCarlo(ngenerators);
  unsigned iters = 10;
  for (unsigned i = 0; i < iters; i++) {
    SurfData* influencers = ab.sampleMonteCarlo(ninfluencers);
    vector<SurfData> closestSets(ngenerators);
    for (unsigned samp = 0; samp < influencers->size(); samp++) {
      unsigned nearest = findClosest(*generators,(*influencers)(samp));
      closestSets[nearest].addPoint((*influencers)[samp]);
    } // for each sample pt
    // Find centroids, update generators
    SurfData* new_generators = new SurfData;
    for (unsigned gen = 0; gen < ngenerators; gen++) {
      if (closestSets[gen].size() != 0) {
        SurfPoint center = computeCentroid(closestSets[gen]);
        double genweight = minalpha + (maxalpha - minalpha)*((double)i/iters);
        new_generators->addPoint(SurfPoint(
          surfpack::weightedAvg((*generators)(gen),center.X(),genweight)));
      } else {
        new_generators->addPoint((*generators)(gen));
      }
    }
    delete generators; generators = new_generators;
    delete influencers;
  } // end iteration
  SurfData result(*generators);
  delete generators;
  return result;
}

VecRbf makeRbfs(const SurfData& generators, const SurfData& radii)
{
  assert(generators.size());
  assert(generators.size() == radii.size());
  vector<RadialBasisFunction> rbfs;
  for (unsigned i = 0; i < generators.size(); i++) {
    rbfs.push_back(RadialBasisFunction(generators(i),radii(i)));
  }
  return rbfs;
}

// Add additional rbfs with broader support to set of candidates
void augment(VecRbf& rbfs)
{
  assert(rbfs.size());
  unsigned toAdd = rbfs.size(); 
  for (unsigned i = 0; i < toAdd; i++) {
    unsigned first = shared_rng()(rbfs.size());
    unsigned second = shared_rng()(rbfs.size());
    //cout << "new basis from " << first << " " << second << endl;
    VecDbl newRadius = rbfs[first].radius;
    if (first == second) { // new function with same center/double radius
      for (unsigned dim = 0; dim < newRadius.size(); dim++) {
        newRadius[dim] *= 2.0;
      }
      rbfs.push_back(RadialBasisFunction(rbfs[first].center,newRadius));
    } else { // new function with avg center, sum of radii
      VecDbl newCenter = surfpack::weightedAvg(rbfs[first].center,rbfs[second].center);
      for (unsigned dim = 0; dim < newRadius.size(); dim++) {
        newRadius[dim] += rbfs[second].radius[dim];
      }
      rbfs.push_back(RadialBasisFunction(newCenter,newRadius));
    }
  }
}

MtxDbl getMatrix(const SurfData& sd, const VecRbf& candidates, VecUns used)
{
  std::sort(used.begin(),used.end());
  MtxDbl A(sd.size(),used.size(),true);
  unsigned nrows = sd.size();
  unsigned ncols = used.size();
  for (unsigned rowa = 0; rowa < nrows; rowa++) {
    for (unsigned cola = 0; cola < ncols; cola++) {
      assert(used[cola] < candidates.size());
      A(rowa,cola) = candidates[used[cola]](sd(rowa));
    }
  }
  return A;
}

VecUns probInclusion(unsigned vec_size, unsigned max_size, double prob)
{
  assert(prob >= 0.0);
  assert(prob <= 1.0);
  assert(vec_size);
  VecUns result;
  for (unsigned i = 0; i < vec_size; i++) {
    if (result.size() >= max_size) break;
    if (shared_rng().randExc() < prob) result.push_back(i);
  }
  return result;
}

VecDbl fullCoeff(unsigned vec_size, const VecDbl& coeffs, VecUns& incl)
{
  VecDbl result(vec_size,0.0);
  for (unsigned i = 0; i < incl.size(); i++) {
    result[incl[i]] = coeffs[i];
  }
  return result;
}

///\todo The exp() function can eat up a lot of time in this method
/// If it becomes a bottleneck, switch to a cached lookup table for
/// the values of exp(x)
//const unsigned granularity = 1000;
//const double maxe = 8.0;
//VecDbl initExps()
//{
//  VecDbl exps;
//  exps.reserve(granularity);
//  for (unsigned i = 0; i < granularity; i++) {
//    exps.push_back(exp(-(double)i/granularity*maxe));
//  }
//  cout << "Initialized " << endl;
//  return exps;
//}
//
//double myexp(const double x)
//{
//  assert(x >= 0.0);
//  static VecDbl exps(initExps());
//  if (x > maxe) return 0.0;
//  return exps[(unsigned)(x/maxe*granularity)];
//}

RadialBasisFunction::RadialBasisFunction(const VecDbl& center_in, const VecDbl& radius_in)
  : center(center_in), radius(radius_in)
{
  assert(!center.empty());
  assert(center.size() == radius.size()); 
}

RadialBasisFunction::RadialBasisFunction(const std::string& center_in, const std::string& radius_in)
  : center(surfpack::toVec<double>(center_in)),
  radius(surfpack::toVec<double>(radius_in))
{
  assert(!center.empty());
  assert(!radius.empty());
  assert(center.size() == radius.size()); 
}

double RadialBasisFunction::operator()(const VecDbl& x) const
{
  assert(x.size() == center.size());
  double sum = 0.0;
  double temp;
  for (unsigned i = 0; i < center.size(); i++) {
    temp = x[i] - center[i];
    sum += temp*temp*radius[i];
  };
  //return myexp(sum);
  return exp(-sum);
}

double RadialBasisFunction::deriv(const VecDbl& x, const VecUns& vars) const
{
  assert(vars.size() == 1);
  assert(!center.empty());
  assert(!radius.empty());
  assert(x.size() == center.size());
  unsigned i = vars[0];
  return -2.0*radius[i]*(x[i]-center[i])*(*this)(x);
}

std::string RadialBasisFunction::asString() const
{
  std::ostringstream os;
  os << "center: ";
  copy(center.begin(),center.end(),std::ostream_iterator<double>(os," "));
  os << " radius: ";
  copy(radius.begin(),radius.end(),std::ostream_iterator<double>(os," "));
  os << std::endl;
  return os.str();
}


RadialBasisFunctionModel::RadialBasisFunctionModel(const VecRbf& rbfs_in, const VecDbl& coeffs_in)
  : SurfpackModel(1), rbfs(rbfs_in),coeffs(coeffs_in)
{
  assert(!rbfs.empty());
  this->ndims = rbfs[0].center.size();
  assert(this->size() != 0);
  assert(rbfs.size() == coeffs.size()); 
}

double RadialBasisFunctionModel::evaluate(const VecDbl& x) const
{
  double sum = 0.0;
  for (unsigned i = 0; i < rbfs.size(); i++) {
    sum += coeffs[i]*rbfs[i](x);
  }
  return sum;
}

/// Currently set up so that operator() must be called immediately before
/// Not good assumption
VecDbl RadialBasisFunctionModel::gradient(const VecDbl& x) const
{
  /// code copied straight from LRM
  assert(!x.empty());
  //assert(coeffs.size() == bs.bases.size());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl result(x.size(),0.0);
  for (unsigned i = 0; i < x.size(); i++) {
    diff_var[0] = i;
    for (unsigned j = 0; j < rbfs.size(); j++) {
      result[i] += coeffs[j]*rbfs[j].deriv(x,diff_var);
    }
  }
  return result;
}

std::string RadialBasisFunctionModel::asString() const
{
  std::ostringstream os;
  unsigned num_bases = rbfs.size();
  unsigned num_vars = ndims;
  os << "-----\n";
  os << "Surfpack Radial Basis Function model\n";
  os << "f(x) = w*phi(x) and phi_k(x) = exp{-r_k*(x-c_k^T).^2}; where\n\n";
  os << "inputs = " << num_vars << "\n";
  os << "bases = " << num_bases << "\n";
 
  os << std::scientific << std::setprecision(16);
  os << "\nw (1 x bases) =\n";
  for(unsigned i=0; i < num_bases; i++) 
    os << std::setw(23) << coeffs[i] << " ";
  os << "\n\nr (bases x inputs) = \n";
  for(unsigned i=0; i < num_bases; i++) {
    for(unsigned j=0; j < num_vars; j++) {
      os << std::setw(23) << rbfs[i].radius[j] << " ";
    }
    os << "\n";
  }
  os << "\nc (bases x inputs) = \n";
  for(unsigned i=0; i < num_bases; i++) {
    for(unsigned j=0; j < num_vars; j++) {
      os << std::setw(23) << rbfs[i].center[j] << " ";
    }
    os << "\n";
  }
  os << "\n-----\n";
  return os.str();
}


typedef std::pair<double,VecUns> RbfBest;
///////////////////////////////////////////////////////////
///	Moving Least Squares Model Factory
///////////////////////////////////////////////////////////

RadialBasisFunctionModelFactory::RadialBasisFunctionModelFactory()
  : SurfpackModelFactory(), nCenters(0), cvtPts(0), maxSubsets(0), 
  minPartition(1)
{

}

RadialBasisFunctionModelFactory::RadialBasisFunctionModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), nCenters(0), cvtPts(0), maxSubsets(0), 
  minPartition(1)
{

}

void RadialBasisFunctionModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["centers"];
  if (strarg != "") nCenters = std::atoi(strarg.c_str());
  strarg = params["cvt_pts"];
  if (strarg != "") cvtPts = std::atoi(strarg.c_str());
  strarg = params["max_subsets"];
  if (strarg != "") maxSubsets = std::atoi(strarg.c_str());
  strarg = params["min_partition"];
  if (strarg != "") minPartition = std::atoi(strarg.c_str());
}

SurfpackModel* RadialBasisFunctionModelFactory::Create(const SurfData& sd)
{
  unsigned max_centers = 100;
  unsigned max_max_subsets = 100;
  if (nCenters == 0) nCenters = min(max_centers,sd.size());
  if (cvtPts == 0) cvtPts = 10*nCenters;
  if (maxSubsets == 0) maxSubsets = min(max_max_subsets,3*nCenters);
  RbfBest bestset(std::numeric_limits<double>::max(),VecUns());
  
  SurfData centers = cvts(AxesBounds::boundingBox(sd),nCenters,cvtPts);
  SurfData radiuses = radii(centers);
  VecDbl b = sd.getResponses();
  VecRbf candidates = makeRbfs(centers,radiuses);
  augment(candidates);
  assert(candidates.size() == 2*nCenters);
  for (unsigned i = 0; i < maxSubsets; i++) {
    VecUns used = probInclusion(candidates.size(),sd.size(),.5);
    MtxDbl A = getMatrix(sd,candidates,used);
    VecDbl x;
    surfpack::linearSystemLeastSquares(A,x,b);
    VecDbl coeffs = fullCoeff(candidates.size(),x,used);
    RadialBasisFunctionModel rbfm(candidates,coeffs);
    StandardFitness sf;
    double fitness = sf(rbfm,sd);
    if (fitness < bestset.first) bestset = RbfBest(fitness,used);
  }
  VecUns used = bestset.second;
  VecRbf final_rbfs;
  VecUns final_used(used.size());
  for (unsigned i = 0; i < used.size(); i++) {
    final_used[i] = i;
    final_rbfs.push_back(candidates[used[i]]);
  }
  // Recompute the coefficients.  If we cached the result, we wouldn't
  // have to do it again.  
  MtxDbl A = getMatrix(sd,final_rbfs,final_used);
  VecDbl x;
  surfpack::linearSystemLeastSquares(A,x,b);
  SurfpackModel* sm = new RadialBasisFunctionModel(final_rbfs, x); 
  StandardFitness sf;
  double fitness = sf(*sm,sd);
  //cout << "Cached fitness: " << bestset.first << " recomputed: " << fitness << endl;
  assert(sm);
  return sm; 
}

