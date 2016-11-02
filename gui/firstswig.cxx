/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "firstswig.h"
#include <iostream>
using namespace std;

FirstClass::FirstClass(int x_) : x(x_) {}

FirstClass::~FirstClass() {}

void FirstClass::printVal(char* msg)
{
  string smsg = string(msg);
  cout << smsg << " " << x++ << endl;
}

double FirstClass::shiftArray(double* vals, int size)
{
  double result = 0;
  for (int i = 0; i < size; i++) 
  {
    result += vals[i];
    vals[i] += x;
  }
  return result;
}

double FirstClass::myevaluate(double* vals, int size)
{
  double result = 0.0;
  result = vals[0] + vals[1];
  result += x;
  return result;
}
