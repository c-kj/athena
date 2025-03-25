// 编译运行：
// icpx -std=c++11 -O3 test_fp.cpp -o test_fp && ./test_fp > test_fp_O3.log

#include <iostream>
#include <cmath>
#include <limits>

#include "../../utils.hpp"

using Real = double;

void test_fp() {
  Real values[] = {
    0.0,
    -0.0,
    std::numeric_limits<Real>::infinity(),
    -std::numeric_limits<Real>::infinity(),
    std::numeric_limits<Real>::quiet_NaN(),
    std::numeric_limits<Real>::denorm_min(),
    std::numeric_limits<Real>::min(),
    std::numeric_limits<Real>::max(),
    1e-300,  // Very small number
    -1e-300, // Very small negative number
    1e300,   // Very large number
    -1e300   // Very large negative number
  };

  for (Real x : values) {
    std::cout << "Testing value: " << x << "\n";

    // 看 fpclassify 的结果
    std::cout << "  -> fpclassify(x) = " << std::fpclassify(x) << "\n";

    // 分类，用标准库函数
    std::cout << "  -> 用标准库函数分类: ";
    if (std::isnan(x)) {
      std::cout << "NaN";
    } else if (std::isinf(x)) {
      std::cout << "inf";
    } else if (x == 0.0) {
      std::cout << "0.0";
    } else if (!std::isnormal(x) && x != 0.0) {
      std::cout << "subnormal";
    } else if (std::isnormal(x)) {
      std::cout << "normal";
    } else {
      std::cout << "CANNOT CLASSIFY THIS VALUE.";
    }
    std::cout << "\n";

    // std::isfinite(x)

    // 分类，用自定义函数
    std::cout << "  -> 用自定义函数分类: ";
    if (fpcheck::is_nan(x)) {
      std::cout << "NaN";
    } else if (fpcheck::is_inf(x)) {
      std::cout << "inf";
    } else if (fpcheck::is_zero(x)) {
      std::cout << "0.0";
    } else if (fpcheck::is_subnormal(x)) {
      std::cout << "subnormal";
    } else if (fpcheck::is_normal(x)) {
      std::cout << "normal";
    } else {
      std::cout << "CANNOT CLASSIFY THIS VALUE.";
    }
    std::cout << "\n";
    
    // Test dangerous operations
    try {
      Real reciprocal = 1.0 / x;
      std::cout << "  -> Reciprocal (1/x): " << reciprocal << "\n";
    } catch (...) {
      std::cout << "  -> Reciprocal (1/x) caused an exception.\n";
    }

    try {
      Real result = x / 0.0;
      std::cout << "  -> Division by zero (x/0): " << result << "\n";
    } catch (...) {
      std::cout << "  -> Division by zero (x/0) caused an exception.\n";
    }

    try{
      Real result = x / x;
      std::cout << "  -> Division by itself (x/x): " << result << "\n";
    } catch (...) {
      std::cout << "  -> Division by itself (x/x) caused an exception.\n";
    }

    std::cout << "\n";
  }
}

int main() {
  test_fp();
  return 0;
}

