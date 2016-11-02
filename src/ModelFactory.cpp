/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "ModelFactory.h"
#include "SurfpackModel.h"
#include "LinearRegressionModel.h"
#include "RadialBasisFunctionModel.h"
#include "DirectANNModel.h"
#include "KrigingModel.h"
#include "MovingLeastSquaresModel.h"
#include "MarsModel.h"

class SurfData;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::istringstream;
using std::string;
using surfpack::shared_rng;

SurfpackModelFactory* ModelFactory::createModelFactory(ParamMap& args)
{
  string type = args["type"];
  SurfpackModelFactory* smf; 
  if (type == "") throw string("Model must declare type");
  else if (type == "polynomial") {
    smf = new LinearRegressionModelFactory(args);
  } else if (type == "mls") {
    smf = new MovingLeastSquaresModelFactory(args);
  } else if (type == "rbf") {
    smf = new RadialBasisFunctionModelFactory(args);
  } else if (type == "kriging") {
    smf = new KrigingModelFactory(args);
  } else if (type == "ann") {
    smf = new DirectANNModelFactory(args);
  } else if (type == "mars") {
    smf = new MarsModelFactory(args);
  } else {
    throw string("Model type requested not recognized");
  }
  // WARNING: the RNG is a static global object, so multiple calls
  // will result in the last taking precedence
  string seedstr = args["seed"];
  if (seedstr != "" ) {
    int seed;
    if (!(istringstream(seedstr) >> seed))
      throw string("Error converting seed to int");
    else
      shared_rng().seed(seed);
  }
  return smf;
}
