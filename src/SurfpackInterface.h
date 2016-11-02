/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_system_headers.h"

class AxesBounds;
class ParsedCommand;
class SurfData;
class SurfpackParser;
class SurfpackModel;

namespace SurfpackInterface
{
  SurfData* LoadData(const std::string& filename);
  SurfData* LoadData(const std::string& filename, unsigned n_predictors,
    unsigned n_responses, unsigned n_cols_to_skip);
  SurfpackModel* LoadModel(const std::string& filename);
  void Save(const SurfpackModel* model, const std::string& filename);
  void Save(const SurfData* data, const std::string& filename);
  SurfpackModel* CreateSurface(const SurfData* sd, ParamMap& args);
  void Evaluate(const SurfpackModel* model, SurfData* sd, 
    const std::string& response_name = "");
  void Evaluate(SurfData* sd, const VecStr test_functions);
  AxesBounds* CreateAxes(const std::string axes);
  SurfData* CreateSample(const AxesBounds* axes, const VecUns grid_points);
  SurfData* CreateSample(const AxesBounds* axes, unsigned n_samples);
  double Fitness(const SurfpackModel*, SurfData* sd, 
    const std::string& metric, unsigned response = 0, unsigned n = 0);
  double Fitness(const SurfpackModel*, const std::string& metric, 
    unsigned response = 0, unsigned n = 0);

  /// whether the interface as compiled supports the passed feature name
  bool HasFeature(const std::string& feature);
};

