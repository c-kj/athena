// Athena++ headers
#include "region.hpp"



bool in_heart(Real a, Real b, Real x, Real y, Real z) {
  x = sqrt(a)*x;
  y = sqrt(a)*y + 0.5;
  return (a-x*x > 0 && pow(x*x, 1./5.) - b*sqrt(a-x*x) < y && y < pow(x*x, 1./5.) + b*sqrt(a-x*x));
}