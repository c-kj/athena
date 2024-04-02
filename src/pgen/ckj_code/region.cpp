
#include <string>     // std::to_string
#include <stdexcept>  // std::invalid_argument

// Athena++ headers
#include "../../athena.hpp"

// 自定义的头文件
#include "region.hpp"



Ball::Ball(const Point &center, const Real radius, const int ndim) : center(center), radius(radius) {
  dim = ndim;
  if (dim == 3) {
    volume = 4./3.*PI*CUBE(radius);
  } else if (dim == 2) {
    volume = PI*SQR(radius);
  } else if (dim == 1) {
    volume = 2*radius;
  } else {
    throw std::invalid_argument("Unknown dimension: " + std::to_string(dim));
  }
}

bool Ball::contains(const Point &point) const {
  if (volume == 0) return false; // 这样当 radius = 0 时，contains() 总是返回 false，避免了除以 0 的问题
  // 根据 dim 是几维，把额外的坐标差都设为 0，只在前 dim 个坐标上计算
  Real dx = dim >= 1 ? point[0] - center[0] : 0;
  Real dy = dim >= 2 ? point[1] - center[1] : 0;
  Real dz = dim >= 3 ? point[2] - center[2] : 0;
  return (dx*dx + dy*dy + dz*dz <= radius*radius); // 用 <= 号，这样方便把临界的 //TODO 这里用 < 还是 <= ？
}



Heart::Heart(const Point &center, const Real size, const Real a, const Real b) : center(center), size(size), a(a), b(b) {
  dim = 2; // 如果要兼容 3D 的话，需要仔细考虑
  volume = PI*a*b * SQR(size);
}

bool Heart::contains(const Point &point) const {
  // 先按 center 和 size 进行平移缩放，然后只需处理标准形状
  Real x = (point[0] - center[0])/size;
  Real y = (point[1] - center[1])/size;
  // 不需要 z，因为是 2D 的

  x = sqrt(a)*x;
  y = sqrt(a)*y + 0.5;
  return (a-x*x > 0 && pow(x*x, 1./5.) - b*sqrt(a-x*x) < y && y < pow(x*x, 1./5.) + b*sqrt(a-x*x));
}