/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "SurfpackInterpreter.h"

using std::vector;
using std::string;

int main(int argc, char** argv)
{
  // If a command line argument is given, us it as the input file
  // The default is to just use standard input
  SurfpackInterpreter si;
  if (argc == 2) {
    string infile(argv[1]);
    si.execute(&infile);
    //infile.close();
  } else {
    si.execute();
  }
  return 0;
}
