/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_c_interface.h"
#include "SurfpackInterface.h"
#include "SurfpackModel.h"

/* Implementation of simplified C interface to some Surfpack library
   functions */

/// one global instance of a SurfpackModel
static SurfpackModel* surfpackCModel = NULL;

// Model load example copied from SurfpackInterpreter.cpp; 
extern "C" 
int surfpack_load_model(const char * const model_filename)
{
  bool success = true;
  try {
    surfpackCModel = SurfpackInterface::LoadModel(std::string(model_filename));
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
double surfpack_eval_model(const double * const eval_pt, unsigned int num_vars)
{
  // evaluate the surfpack model, using a std::vector
  std::vector<double> eval_vec(eval_pt, eval_pt + num_vars);
  return surfpackCModel->operator()(eval_vec);
}


extern "C"
void surfpack_free_model()
{
  delete surfpackCModel;
}
