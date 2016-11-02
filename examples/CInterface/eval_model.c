/* Example of loading a Surfpack model from a surface file (.sps or
   .bsps) and evaluating it at a number of points. */

#include <stdlib.h>
#include <stdio.h>
#include "surfpack_c_interface.h"


int main()
{
  const char* const model_name = "sp_gp_model.sps";
  const int num_vars = 2;

  if (!surfpack_load_model(model_name)) {

    // evaluate the model at the point, with a perturbation
    double base_pt[2] = {0.97, 0.94};
    double eval_pt[num_vars];

    // perform one million trials
    int num_trials = 5;
    int i;
    for (i=0; i<num_trials; ++i) {

      int j;
      // add +/-1% random perturbation to all components
      for (j=0; j<num_vars; ++j) {
	eval_pt[j] = base_pt[j] + 
	  -0.01 + 0.02 * (double) rand() / (double) RAND_MAX;
      }

      // evaluate the surfpack model
      double resp = surfpack_eval_model(eval_pt, num_vars);
      printf("Response is: %g\n", resp);

    }

    surfpack_free_model();
    return 0;

  }
  else
    printf("Error loading model\n");

  return 1;
}
