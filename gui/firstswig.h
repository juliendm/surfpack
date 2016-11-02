/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

class FirstClass
{
public:
  FirstClass(int x_);
  ~FirstClass();
  void printVal(char* msg);
  double shiftArray(double* vals, int size);
  double myevaluate(double* vals, int size);
private:
  int x;
};
