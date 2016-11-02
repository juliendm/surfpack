/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef FLEX_WRAPPER_H
#define FLEX_WRAPPER_H

#include "surfpack_system_headers.h"

class FlexWrapper
{
public:
  FlexWrapper();
  ~FlexWrapper();
  void setParseStreams(const std::string* input_string, const std::string* output_string);
  int nextToken();
  const char* currentToken();
private:
  FILE* infile;
  FILE* outfile;
};

#endif
