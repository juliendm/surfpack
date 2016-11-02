/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_C_INTERFACE_H
#define SURFPACK_C_INTERFACE_H

/* Simplified C interface to Surfpack libraries to load and evaluate a
   previously saved model.  Uses a single static instance of a
   Surfpack model, so cannot be used for concurrent models in a single
   program.

   Owner: Brian M. Adams, briadam@sandia.gov

   TODO: Give the client to struct for the model, so the client can
   manage multiple concurrent model instances from C.
*/


/* Load a surfpack model from the specified .sps or .bsps file.
   Returns 0 on success, 1 on failure. */
#ifdef __cplusplus
extern "C" 
#endif
int surfpack_load_model(const char * const model_filename);

/* Evaluate the loaded model at the provided eval_pt, with length
   num_vars.  Returns the surrogatem model value at eval_pt. */
#ifdef __cplusplus
extern "C"
#endif
double surfpack_eval_model(const double * const eval_pt, unsigned int num_vars);

/* Free the surfpack model from memory. */
#ifdef __cplusplus
extern "C"
#endif
void surfpack_free_model();


#endif  /* SURFPACK_C_INTERFACE_H */
