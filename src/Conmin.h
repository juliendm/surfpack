/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __CONMIN_H__ 
#define __CONMIN_H__ 

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"


class Conmin {
public:
  Conmin(unsigned ndv_in);
  void bounds(const VecDbl& lower_bounds, const VecDbl& upper_bounds);
  virtual void optimize(VecDbl& x, double& final_val, unsigned max_iter) = 0;
  virtual double objective(const VecDbl& x) = 0;
  virtual VecDbl gradient(const VecDbl& x) = 0;
  virtual ~Conmin();
protected:
  VecDbl upperBounds;
  VecDbl lowerBounds;
  int NSIDE;
  unsigned ndv;
  
};

#endif
