/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

//#include "StdAfx.h"
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include "unittests.h"
#include "SurfaceTest.h"
#include "SurfDataTest.h"
#include "SurfScalerTest.h"
//#include "PolynomialSurfaceTest.h"
#include "SurfPointTest.h"
#include "SurfpackCommonTest.h"
#include "LinearRegressionModelTest.h"
#include "MovingLeastSquaresTest.h"
#include "RadialBasisFunctionTest.h"
#include "SurfpackModelTest.h"
#include "KrigingModelTest.h"
#include "ModelScalerTest.h"
#include "ModelFactoryTest.h"
#include <string>

int main(int argc, char* argv[])
{
  // Get the top level suite from the registry
  //CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
  CppUnit::Test *surfpoint_suite = SurfPointTest::suite(); 
  CppUnit::Test *surf_data_suite = SurfDataTest::suite(); 
  CppUnit::Test *surf_scaler_suite = SurfScalerTest::suite(); 
  CppUnit::Test *surface_suite = SurfaceTest::suite(); 
  //CppUnit::Test *polynomial_suite = PolynomialSurfaceTest::suite(); 
  CppUnit::Test *surfpack_common_suite = SurfpackCommonTest::suite(); 
  CppUnit::Test *linear_regression_model_suite = LinearRegressionModelTest::suite(); 
  CppUnit::Test *surfpack_model_suite = SurfpackModelTest::suite(); 
  CppUnit::Test *model_scaler_suite = ModelScalerTest::suite(); 
  CppUnit::Test *moving_least_squares_suite = MovingLeastSquaresTest::suite(); 
  CppUnit::Test *radial_basis_function_suite = RadialBasisFunctionTest::suite(); 
  CppUnit::Test *kriging_model_test_suite = KrigingModelTest::suite(); 
  CppUnit::Test *model_factory_test_suite = ModelFactoryTest::suite(); 

  // Adds the test to the list of test to run
  CppUnit::TextUi::TestRunner runner;
  //runner.addTest( surfpoint_suite );
  //runner.addTest( surf_data_suite );
  //runner.addTest( surf_scaler_suite );
  //runner.addTest( surface_suite );
  //runner.addTest( polynomial_suite );
  //runner.addTest( surfpack_common_suite );
  //runner.addTest( linear_regression_model_suite );
  //runner.addTest( surfpack_model_suite );
  //runner.addTest( moving_least_squares_suite );
  //runner.addTest( radial_basis_function_suite );
  //runner.addTest( model_scaler_suite );
  //runner.addTest( kriging_model_test_suite );
  runner.addTest( model_factory_test_suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
                                                       std::cerr ) );
  //runner.setOutputter( new CppUnit::XmlOutputter( &runner.result(),
  //                                                     std::cerr ) );
  // Run the test.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of tests failed.
  return wasSucessful ? 0 : 1;
}
