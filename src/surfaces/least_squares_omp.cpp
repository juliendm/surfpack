/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "least_squares_omp.h"
#ifdef HAVE_PECOS
#include "pecos_data_types.hpp"  // for basic scalar and matrix data types
#include "CrossValidation.hpp"   // for ValidationIterator
#include "LinearSolver.hpp"      // for CompressedSensingOptions/Tool
#endif

namespace surfpack {

/** Input args can't be const due to potential forward to
    linearSystemLeastSquares() */
void leastSquaresOMP(MtxDbl& A_in, VecDbl& b_in, int random_seed, VecDbl& x_out)
{

#ifdef HAVE_PECOS
  int num_data = b_in.size();
  int num_coeffs = A_in.getNCols();

  if (A_in.getNRows() != num_data) {
    std::cerr << "\nError (Surfpack OMP): linear system size mismatch" 
	      << std::endl;
    throw std::runtime_error("Surfpack matrix size mismatch");
  }

  // copy the data to Teuchos for use with pecos
  Pecos::RealMatrix A(A_in.getNRows(), A_in.getNCols());
  Pecos::RealVector b(b_in.size());
  for(size_t r=0; r<A_in.getNRows(); ++r) {
    for(size_t c=0; c<A_in.getNCols(); ++c)
      A(r,c) = A_in(r,c);
    b[r] = b_in[r]; 
  }

  try {
 
    Pecos::MultipleSolutionLinearModelCrossValidationIterator cv_iterator;

    // seed for CV splitting
    cv_iterator.set_seed(random_seed);

    // TODO: NaN / Inf
    // cv_iterator.set_fault_data( faultInfo,
    // 			      surrData.failed_response_data() );

    Pecos::CompressedSensingOptions CSOpts;
    CSOpts.solver = Pecos::ORTHOG_MATCH_PURSUIT;

    Pecos::CompressedSensingTool CSTool;
    CSTool.set_linear_solver(CSOpts);
  
    Pecos::LinearSolver_ptr linear_solver = CSTool.get_linear_solver();

    cv_iterator.set_solver(linear_solver);

    // perform at most 10 cross validation folds (surfpack models may
    // be trained on limited data)
    int num_folds = std::min(num_data, 10);

    cv_iterator.set_max_num_unique_tolerances(100);
    cv_iterator.set_num_folds(num_folds);
    cv_iterator.set_num_points(num_data);
    // if ( use_gradients )
    //   cv_iterator.set_num_equations_per_point( sharedDataRep->numVars + 1 );
    // else 
    cv_iterator.set_num_equations_per_point(1);

    Pecos::Real score = cv_iterator.run_cross_validation(A, b);
    Pecos::Real best_tolerance = cv_iterator.get_best_residual_tolerance();

    Pecos::RealMatrix solutions, metrics;
    linear_solver->set_residual_tolerance(best_tolerance);
    linear_solver->solve(A, b, solutions, metrics);

    Pecos::Real* dense_coeffs = solutions[solutions.numCols()-1];
    x_out.resize(num_coeffs);
    for(size_t r=0; r<num_coeffs; ++r)
      x_out[r] = dense_coeffs[r];

  }
  catch (const std::exception& e) {
    std::cerr << "\nError Surfpack exception in OMP:\n" << e.what() 
	      << std::endl;
    throw;
  }
#else
  surfpack::linearSystemLeastSquares(A_in,x_out,b_in);
#endif

}

}  // namespace surfpack
