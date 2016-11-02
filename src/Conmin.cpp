#include "Conmin.h"


using std::cout;
using std::endl;
using std::vector;

Conmin::Conmin(unsigned ndv_in)
  : NSIDE(0), ndv(ndv_in)
{
  cout << "ndv: " << ndv << endl;
  assert(ndv > 0);
  const unsigned maxdims = 40;
  assert(ndv <= maxdims);
}

void Conmin::bounds(const VecDbl& lower_bounds, const VecDbl& upper_bounds)
{
  assert(upper_bounds.size() == lower_bounds.size());
  upperBounds = upper_bounds;
  lowerBounds = lower_bounds;
  NSIDE=1; // I believe this specifies that the side constraints are l/u bounds
}

Conmin::~Conmin()
{

}

