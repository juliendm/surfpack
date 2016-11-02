/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#ifndef MODEL_FACTORY_TEST_H 
#define MODEL_FACTORY_TEST_H 

#include <cppunit/extensions/HelperMacros.h>

class ModelFactoryTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ModelFactoryTest );
CPPUNIT_TEST( simpleTest );
CPPUNIT_TEST( argsTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
void simpleTest();
void argsTest();
};

#endif
