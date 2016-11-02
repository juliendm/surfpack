/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_system_headers.h"
#include "SurfpackInterface.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackParserArgs.h"
#include "SurfData.h"
#include "ModelFactory.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;
using std::ostringstream;

#include "SurfpackModel.h"
#include "ModelFitness.h"
using SurfpackInterface::CreateAxes;
using SurfpackInterface::CreateSurface;
using SurfpackInterface::LoadData;
using SurfpackInterface::LoadModel;
using SurfpackInterface::Evaluate;

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterface namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////


SurfData* SurfpackInterface::LoadData(const std::string& filename)
{
  SurfData* sd = new SurfData(filename);
  assert(sd);
  return sd;
}

SurfData* SurfpackInterface::LoadData(const std::string& filename, unsigned n_predictors, unsigned n_responses, unsigned n_cols_to_skip)
{
  SurfData* sd = new SurfData(filename,n_predictors,n_responses,n_cols_to_skip);
  assert(sd);
  return sd;
}

SurfpackModel* SurfpackInterface::LoadModel(const std::string& filename)
{
  SurfpackModel* model = NULL;

  // TODO: clean up where files get opened / closed
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  bool binary = surfpack::isBinaryModelFilename(filename);
  
  std::ifstream model_ifstream(filename.c_str(), (binary ? std::ios::in|std::ios::binary : std::ios::in));
  if (!model_ifstream.good())
    throw std::string("Failure opening model file for load."); 

  if (binary) {
    boost::archive::binary_iarchive input_archive(model_ifstream);
    input_archive >> model; 
    std::cout << "Model loaded from binary file '" << filename << "'." 
	      << std::endl;
  }
  else {
    boost::archive::text_iarchive input_archive(model_ifstream);
    input_archive >> model; 
    std::cout << "Model loaded from text file '" << filename << "'." 
	      << std::endl;
  }

#else
  throw string("surface load requires compilation with Boost serialization.");
#endif

  return model;
}

void SurfpackInterface::Save(const SurfpackModel* model, const std::string& filename)
{
  // TODO: consider where files are opened/managed
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  bool binary = surfpack::isBinaryModelFilename(filename);

  std::ofstream model_ofstream(filename.c_str(), (binary ? std::ios::out|std::ios::binary : std::ios::out));  
  if (!model_ofstream.good())
    throw std::string("Failure opening model file for save."); 

  if (binary) {
    boost::archive::binary_oarchive output_archive(model_ofstream);
    output_archive << model;
    std::cout << "Model saved to binary file '" << filename << "'." 
	      << std::endl;
  }
  else {
    boost::archive::text_oarchive output_archive(model_ofstream);
    output_archive << model;
    std::cout << "Model saved to text file '" << filename << "'." << std::endl;
  }
#else
  throw 
    string("surface save requires compilation with Boost serialization.");
#endif
}




void SurfpackInterface::Save(const SurfData* data, const std::string& filename)
{
  data->write(filename);
}

SurfpackModel* SurfpackInterface::CreateSurface(const SurfData* sd, ParamMap& args)
{
  assert(sd);
  SurfpackModel* model = 0;
  // Surface* model = ModelFactory::Create(ParamMap)
  return model;
}

void SurfpackInterface::Evaluate(const SurfpackModel* model, SurfData* sd, 
  const std::string& response_name)
{
  assert(model);
  assert(sd);
  VecDbl responses = (*model)(*sd);
  sd->addResponse(responses, response_name);
}

void SurfpackInterface::Evaluate(SurfData* sd, const VecStr test_functions)
{
  assert(sd);
  for (VecStr::const_iterator itr = test_functions.begin();
    itr != test_functions.end(); ++itr) {
    VecDbl results(sd->size());
    for (unsigned i = 0; i < results.size(); i++) {
      results[i] = surfpack::testFunction(*itr,(*sd)(i));
    }
    sd->addResponse(results,*itr);
  } 
}

AxesBounds* SurfpackInterface::CreateAxes(const std::string axes)
{
  return new AxesBounds(axes);
}

SurfData* SurfpackInterface::CreateSample(const AxesBounds* axes, const VecUns grid_points)
{
  return axes->sampleGrid(grid_points);  
}

SurfData* SurfpackInterface::CreateSample(const AxesBounds* axes, unsigned n_samples)
{
  return axes->sampleMonteCarlo(n_samples);
}

double SurfpackInterface::Fitness(const SurfpackModel* model, SurfData* sd, 
const std::string& metric, unsigned response, unsigned n)
{
  assert(model);
  assert(sd);
  sd->setDefaultIndex(response);
  ModelFitness* mf = ModelFitness::Create(metric,n);
  double result = (*mf)(*model,*sd);
  delete mf;
  return result;
}

/// Doxygen comment
double SurfpackInterface::Fitness(const SurfpackModel*, const std::string& metric, 
unsigned response, unsigned n)
{
  throw string("Must pass data set to compute metric");
  return 0.0;
}

bool SurfpackInterface::HasFeature(const std::string& feature)
{
  if (feature == "model_save" || feature == "model_load") {
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    return true;
#else
    return false;
#endif
  }

  return false;
}

