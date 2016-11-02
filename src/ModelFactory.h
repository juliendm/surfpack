/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Model Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MODEL_FACTORY_H__
#define __MODEL_FACTORY_H__

#include "surfpack_system_headers.h"

class SurfpackModelFactory;
class SurfpackModel;
/// The createModel methods are intended to be a sort of virtual constructor.
/// When new Model sub-classes are added, the changes can be made here 
/// without having to touch the Model class itself.
/// \todo Expand the ModelFactory namespace into a singleton class.  Add 
/// pairs of strings and function pointers to an STL map so that the 
/// createModel methods just do a lookup in the map instead of the clunky
/// if...else construct.  Priority: low.
namespace ModelFactory {
  SurfpackModelFactory* createModelFactory(ParamMap& args);
}
#endif
