/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_LSQ_OMP
#define SURFPACK_LSQ_OMP

#include "surfpack.h"

namespace surfpack {

/// Least squares solver performing basis selection using Pecos
/// implementation of cross validation-guided orthogonal matching pursuit.
/// If pecos not available, reverts to standard surfpack LSQ solver
void leastSquaresOMP(MtxDbl& A_in, VecDbl& b_in, int random_seed, VecDbl& x_out);

}
#endif
