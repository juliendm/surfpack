#include "surfpack_system_headers.h"
#include "ModelScaler.h"
#include "SurfData.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(ModelScaler)
BOOST_CLASS_EXPORT(NonScaler)
BOOST_CLASS_EXPORT(NormalizingScaler)
#endif


const VecDbl& NonScaler::scale(const VecDbl& unscaled_x) const 
{ 
  return unscaled_x;
}

double NonScaler::descale(double scaled_response) const 
{
  return scaled_response;
}

double NonScaler::scaleResponse(double unscaled_response) const 
{
  return unscaled_response; 
}

ModelScaler* NonScaler::Create(const SurfData& data)
{
  return new NonScaler(); 
}

std::string NonScaler::asString()
{
  return string("No scaling");
}

ModelScaler* NonScaler::clone() const
{
  return new NonScaler(*this);
}

const VecDbl& NormalizingScaler::scale(const VecDbl& unscaled_x) const
{
  //cout << "NormalizingScaler::scale" << endl;
  if(unscaled_x.size() != scalers.size()) {
    std::cout << "unscaled_x.size=" << unscaled_x.size() <<
      " scalers.size=" << scalers.size() << std::endl;
    assert(unscaled_x.size() == scalers.size());
  }
  assert(this->result.size() == scalers.size());
  for (unsigned i = 0; i < scalers.size(); i++) {
    this->result[i] = (unscaled_x[i] - scalers[i].offset)/scalers[i].scaleFactor;
  }
  return this->result;
}

double NormalizingScaler::descale(double scaled_response) const
{
  return scaled_response*descaler.scaleFactor + descaler.offset;
}

double NormalizingScaler::scaleResponse(double unscaled_response) const 
{
  return (unscaled_response - descaler.offset) / descaler.scaleFactor;
}


VecDbl NormalizingScaler::getScalerOffsets() const
{
  VecDbl offsets(scalers.size());
  for(unsigned i = 0; i < scalers.size(); i++) {
    offsets[i] = scalers[i].offset;
  }
  return offsets;
}

VecDbl NormalizingScaler::getScalerScaleFactors() const
{
  VecDbl scaleFactors(scalers.size());
  for(unsigned i = 0; i < scalers.size(); i++) {
    scaleFactors[i] = scalers[i].scaleFactor;
  }
  return scaleFactors;
}

double NormalizingScaler::getDescalerOffset() const
{
  return descaler.offset;
}

double NormalizingScaler::getDescalerScaleFactor() const
{
  return descaler.scaleFactor;
}

ModelScaler* NormalizingScaler::Create(const SurfData& data)
{
  vector<NormalizingScaler::Scaler> scalers(data.xSize());
  for (unsigned i = 0; i < data.xSize(); i++) {
    VecDbl predictor = data.getPredictor(i);
    scalers[i].offset = *(std::min_element(predictor.begin(),predictor.end()));
    scalers[i].scaleFactor = 
      *(std::max_element(predictor.begin(),predictor.end())) - scalers[i].offset;
  }
  NormalizingScaler::Scaler descaler;
  VecDbl response = data.getResponses();
  descaler.offset = *(std::min_element(response.begin(),response.end()));
  descaler.scaleFactor = 
    *(std::max_element(response.begin(),response.end())) - descaler.offset;
  return new NormalizingScaler(scalers,descaler);
}

ModelScaler* NormalizingScaler::Create(const SurfData& data,
				       double norm_factor)
{
  double min_elem, max_elem;
  assert(norm_factor >= 0.0);
  vector<NormalizingScaler::Scaler> scalers(data.xSize());
  for (unsigned i = 0; i < data.xSize(); i++) {
    VecDbl predictor = data.getPredictor(i);
    min_elem = *(std::min_element(predictor.begin(),predictor.end()));
    max_elem = *(std::max_element(predictor.begin(),predictor.end()));
    scalers[i].offset = (max_elem + min_elem)/2.0;
    scalers[i].scaleFactor =  (max_elem - min_elem)/2.0/norm_factor;
  }
  NormalizingScaler::Scaler descaler;
  VecDbl response = data.getResponses();
  min_elem = *(std::min_element(response.begin(),response.end()));
  max_elem = *(std::max_element(response.begin(),response.end()));
  descaler.offset = (max_elem + min_elem)/2.0;
  descaler.scaleFactor = (max_elem - min_elem)/2.0/norm_factor;
  return new NormalizingScaler(scalers,descaler);
}

ModelScaler* NormalizingScaler::clone() const
{
  return new NormalizingScaler(*this);
}

std::string NormalizingScaler::asString()
{
  std::ostringstream os;
  for (unsigned i = 0; i < scalers.size(); i++) {
    os << "offset: " << scalers[i].offset 
       << " scaleFactor: " << scalers[i].scaleFactor
       << "\n";
  }
  os << "descaler offset: " << descaler.offset << " scaleFactor: " << descaler.scaleFactor << endl;
  return os.str();
}

/// ScaledSurfData
ScaledSurfData::ScaledSurfData(const ModelScaler& ms_in, const SurfData& sd_in)
  : ms(ms_in), sd(sd_in)
{

}

VecDbl ScaledSurfData::getResponses() const
{
  VecDbl responses = sd.getResponses();
  for (VecDbl::iterator it = responses.begin(); it != responses.end(); ++it) {
    *it = ms.scaleResponse(*it);
  }
  return responses;
}

double ScaledSurfData::getResponse(unsigned index) const
{
  //checkRangeNumPoints(header, index);
  double unscaled = sd.getResponse(index);
  return ms.scaleResponse(unscaled);
}

unsigned ScaledSurfData::size() const
{
  return sd.size();
}

unsigned ScaledSurfData::xSize() const
{
  return sd.xSize();
}

double ScaledSurfData::operator()(unsigned pt, unsigned dim) const
{
  assert(pt < sd.size());
  assert(dim < sd.xSize());
  //return points[mapping[pt]]->X()[dim];
  const VecDbl& unscaled_pt = sd[pt].X();
  const VecDbl& scaled_pt = ms.scale(unscaled_pt);
  return scaled_pt[dim];
}

const VecDbl& ScaledSurfData::operator()(unsigned pt) const
{
  assert(pt < sd.size());
  const VecDbl& unscaled_pt = sd[pt].X();
  return ms.scale(unscaled_pt);
  //copy(scaled_pt.begin(),scaled_pt.end(),std::ostream_iterator<double>(cout," "));
}

VecVecDbl ScaledSurfData::asVecVecDbl(const ScaledSurfData& data)
{
  VecVecDbl result(data.size());
  for (unsigned i = 0; i < data.size(); i++) {
    result[i].resize(data.xSize());
    for (unsigned j = 0; j < data.xSize(); j++) {
      result[i][j] = data(i,j);
    }
  }
  return result;
}
