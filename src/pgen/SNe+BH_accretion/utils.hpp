// 自定义的工具函数

#ifndef CKJ_CODE_UTILS_HPP_
#define CKJ_CODE_UTILS_HPP_

// 自定义函数所需的标准库导入
#include <fstream>  // std::ofstream
#include <map>
#include <chrono>
#include <vector>

// check_nan 函数所需的标准库导入
#include <cstdint>     // for uint32_t, uint64_t
#include <cstring>     // for std::memcpy
#include <type_traits> // for std::is_floating_point

// Athena++ headers
#include "../../athena.hpp"
#include "../../parameter_input.hpp"





void check_and_create_directory(const std::string &path);

void ensure_parent_directory_exists(const std::string &path);

// 泛化的 writeMapToFile 函数，适用于装有各种类型的 map
template <typename K, typename V>
void writeMapToFile(const std::map<K, V>& map, const std::string& filename) {
  ensure_parent_directory_exists(filename);  // 确保父目录存在

  std::ofstream file{filename};
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << '\n';
    return;
  }

  for (const auto& pair : map) {
    file << pair.first << ": " << pair.second << '\n';
  }

  file.close();
}


std::tuple<std::vector<std::array<Real, 3>>, std::vector<int>> get_AMR_points_and_levels(ParameterInput *pin);

// 从字符串中读取数字列表
// 能解析 "1"、"1,2,3"、"1., 2.0, .3,"，但不能 "1,2,3,,"
std::vector<Real> read_vector(std::string str, char delimiter=',');

// 从字符串中读取数字 array，用于 Point、Vector 等
template <size_t N>
std::array<Real, N> read_array(const std::string& str, char delimiter=',') {
  std::array<Real, N> arr;
  std::vector<Real> vec = read_vector(str, delimiter);
  if (vec.size() != N) {
    throw std::invalid_argument("The input string should contain exactly " + std::to_string(N) + " numbers separated by " + delimiter);
  }
  for (size_t i = 0; i < N; i++) {
    arr[i] = vec[i];
  }
  return arr;
}


std::string format_duration(const std::chrono::duration<double>& duration);


// 使用位运算对浮点数 (fp) 进行各类检查、分类
namespace fpcheck {

  // 使用位运算检测浮点数类别，返回标准宏值。保证不受编译器优化影响。
  // 适用于各种系统、平台、编译器（只要符合 IEEE 754 标准）
  // 返回值为标准库中定义的宏值：FP_NAN, FP_INFINITE, FP_ZERO, FP_SUBNORMAL, FP_NORMAL
  template<typename T>
  int fpclassify_bits(T value) {
    static_assert(std::is_floating_point<T>::value, "fpclassify_bits requires floating point type");  // 检查 T 是否为浮点类型
    
    // 根据浮点类型选择对应的整数类型，用于存储位表示
    constexpr bool is_float = sizeof(T) == sizeof(uint32_t);
    using UInt = typename std::conditional<is_float, uint32_t, uint64_t>::type;
    
    // 指数和尾数的掩码
    constexpr UInt exponent_mask = is_float ? 0x7F800000u : 0x7FF0000000000000u;
    constexpr UInt mantissa_mask = is_float ? 0x007FFFFFu : 0x000FFFFFFFFFFFFFu;
    
    // 获取浮点数的位表示
    UInt bits;                                           // 用于存储 value 的位表示
    std::memcpy(&bits, &value, sizeof(T));  // 将 value 的位表示拷贝到 bits 中
    
    // 提取指数和尾数部分
    UInt exponent = bits & exponent_mask;
    UInt mantissa = bits & mantissa_mask;
    
    // 分类判断，使用标准库中定义的宏
    if (exponent == exponent_mask) {  // 指数位全1
      return (mantissa != 0) ? FP_NAN : FP_INFINITE;
    } else if (exponent == 0) {  // 指数位全0
      return (mantissa != 0) ? FP_SUBNORMAL : FP_ZERO;
    } else {  // 其他情况为正规数
      return FP_NORMAL;  // 指数位不为全 1（排除无穷大和 NaN），也不为全 0（排除次正规数和零）。
    }
  }

  //* 以下命名有意加了一个 _ 以与标准库中的区分

  template<typename T>
  bool is_nan(T value) { return fpclassify_bits(value) == FP_NAN; }

  template<typename T>
  bool is_inf(T value) { return fpclassify_bits(value) == FP_INFINITE; }

  template<typename T>
  bool is_zero(T value) { return fpclassify_bits(value) == FP_ZERO; }

  template<typename T>
  bool is_subnormal(T value) { return fpclassify_bits(value) == FP_SUBNORMAL; }

  template<typename T>
  bool is_normal(T value) { return fpclassify_bits(value) == FP_NORMAL; }

  template<typename T>
  bool is_finite(T value) {
    int classification = fpclassify_bits(value);
    return classification != FP_INFINITE and classification != FP_NAN;
  }


} // namespace fpcheck

#endif // CKJ_CODE_UTILS_HPP_