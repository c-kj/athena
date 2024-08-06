#ifndef REGION_HPP_
#define REGION_HPP_

// C++ headers
#include <array>
#include <iostream>

// Athena++ headers
#include "../../athena.hpp"


using Point = std::array<Real, 3>; // 目前把所有情况都当 3D 处理。对于 2D，第三个坐标可以任意。
using Vector = std::array<Real, 3>;


// 重载 << 操作符。实现在 region.cpp 中，从而可以被其他文件多处 include
std::ostream& operator<<(std::ostream& os, const Point& point);

// 矢量的叉乘
inline Vector CrossProduct(const Vector &a, const Vector &b) {
  return {
    a[1]*b[2] - a[2]*b[1], 
    a[2]*b[0] - a[0]*b[2], 
    a[0]*b[1] - a[1]*b[0]
  };
}


struct Region {
  int dim;
  Real volume;

  virtual bool contains(const Point &point) const = 0; // 纯虚函数，必须在派生类中实现
  virtual ~Region() = default; // 虚析构函数，这样当智能指针销毁时，会调用派生类的析构函数，而非基类的（默认空）析构函数
};


struct Ball : public Region {
  Point center;
  Real radius;

  Ball(const Point &center, const Real radius, const int dim);

  bool contains(const Point &point) const override;
};

struct Heart : public Region {
  Point center;
  Real size;
  Real a, b;

  Heart(const Point &center, const Real size, const Real a=3.3, const Real b=0.75);

  bool contains(const Point &point) const override;
};



#endif // REGION_HPP_