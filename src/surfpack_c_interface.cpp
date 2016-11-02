/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_c_interface.h"
#include "SurfpackInterface.h"
#include "SurfData.h"
#include "SurfpackModel.h"
#include "ModelFactory.h"

/* Implementation of simplified C interface to some Surfpack library
   functions */

/// one global instance of a SurfpackModel
//typedef std::pair<std::string, SurfData*> SurfDataSymbol;
typedef std::map<std::string, SurfData*> SurfDataMap;
//typedef std::pair<std::string, SurfpackModel*> SurfpackModelSymbol;
typedef std::map<std::string, SurfpackModel*> SurfpackModelMap;

static SurfDataMap dataVars;
static SurfpackModelMap modelVars;

// Model load example copied from SurfpackInterpreter.cpp; 
extern "C" 
int surfpack_load_data(const char * const name, const char * const data_filename, unsigned int n_predictors, unsigned int n_responses, unsigned int n_cols_to_skip)
{
  bool success = true;
  try {
    SurfData* data = SurfpackInterface::LoadData(std::string(data_filename),n_predictors,n_responses,n_cols_to_skip);
    //dataVars.insert(SurfDataSymbol(std::string(name), data));
    if (dataVars.find(std::string(name)) != dataVars.end()) { delete dataVars[std::string(name)]; }
    dataVars[std::string(name)] = data;
  }
  catch (const std::exception& e) {
    std::cerr << "Error loading surfpack data! Exception:\n" << e.what() 
        << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error loading surfpack data! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error loading surfpack data! Unknown error." << std::endl;
    success = false;
  }
  return (success ? 0 : 1);
}

extern "C" 
int surfpack_save_data(const char * const name, const char * const data_filename)
{
  bool success = true;
  try {
    SurfData* data = dataVars[std::string(name)];
    SurfpackInterface::Save(data, std::string(data_filename));
  }
  catch (const std::exception& e) {
    std::cerr << "Error saving surfpack data! Exception:\n" << e.what() 
        << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error saving surfpack data! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error saving surfpack data! Unknown error." << std::endl;
    success = false;
  }
  return (success ? 0 : 1);
}

extern "C" 
void surfpack_add_data(const char * const name, const double * const build_pt, unsigned int num_vars, const double f)
{
  std::vector<double> build_vec(build_pt, build_pt + num_vars);
  if (dataVars.find(std::string(name)) == dataVars.end()) { dataVars[std::string(name)] = new SurfData(); }
  SurfData* data = dataVars[std::string(name)];
  data->addPoint(SurfPoint(build_vec,f));
}

extern "C" 
int surfpack_load_model(const char * const name, const char * const model_filename)
{
  bool success = true;
  try {
    SurfpackModel* model = SurfpackInterface::LoadModel(std::string(model_filename));
    //modelVars.insert(SurfpackModelSymbol(std::string(name), model));
    if (modelVars.find(std::string(name)) != modelVars.end()) { delete modelVars[std::string(name)]; }
    modelVars[std::string(name)] = model;
  }
  catch (const std::exception& e) {
    std::cerr << "Error loading surfpack model! Exception:\n" << e.what() 
	      << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error loading surfpack model! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error loading surfpack model! Unknown error." << std::endl;
    success = false;
  }
  return (success ? 0 : 1);
}

extern "C" 
int surfpack_save_model(const char * const name, const char * const model_filename)
{
  bool success = true;
  try {
    SurfpackModel* model = modelVars[std::string(name)];
    SurfpackInterface::Save(model, std::string(model_filename));
  }
  catch (const std::exception& e) {
    std::cerr << "Error saving surfpack model! Exception:\n" << e.what() 
        << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error saving surfpack model! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error saving surfpack model! Unknown error." << std::endl;
    success = false;
  }
  return (success ? 0 : 1);
}

extern "C" 
int surfpack_build_model(const char * const name, const char * const type)
{
  bool success = true;
  try {
    SurfData* data = dataVars[std::string(name)];
    std::map<std::string, std::string> args;
    args["type"] = std::string(type);
    SurfpackModelFactory* smf = ModelFactory::createModelFactory(args);
    SurfpackModel* model = smf->Build(*data);
    delete smf;
    if (modelVars.find(std::string(name)) != modelVars.end()) { delete modelVars[std::string(name)]; }
    modelVars[std::string(name)] = model;

  }
  catch (const std::exception& e) {
    std::cerr << "Error building surfpack model! Exception:\n" << e.what() 
        << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error building surfpack model! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error building surfpack model! Unknown error." << std::endl;
    success = false;
  }
  return (success ? 0 : 1);
}


extern "C"
double surfpack_eval_model(const char * const name, const double * const eval_pt, unsigned int num_vars)
{
  std::vector<double> eval_vec(eval_pt, eval_pt + num_vars);
  SurfpackModel* model = modelVars[std::string(name)];
  return model->operator()(eval_vec);
}

extern "C"
double surfpack_variance_model(const char * const name, const double * const eval_pt, unsigned int num_vars)
{
  std::vector<double> eval_vec(eval_pt, eval_pt + num_vars);
  SurfpackModel* model = modelVars[std::string(name)];
  return model->variance(eval_vec);
}

extern "C"
void surfpack_free_model()
{
  for (SurfpackModelMap::iterator iter = modelVars.begin();
       iter != modelVars.end();
       ++iter) {
    delete iter->second; 
  }
}
