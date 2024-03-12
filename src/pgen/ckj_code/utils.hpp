// 自定义的工具函数

#ifndef CKJ_CODE_UTILS_HPP_
#define CKJ_CODE_UTILS_HPP_

// 自定义函数所需的标准库导入
#include <fstream>  // std::ofstream

// Athena++ headers
#include "../../athena.hpp"
#include "../../parameter_input.hpp"





void printf_and_save_to_stream(std::ostream &stream, const char *format, ...);

void check_and_create_directory(const std::string &path);

void ensure_parent_directory_exists(const std::string &path);

std::tuple<std::vector<std::array<Real, 3>>, std::vector<int>> get_AMR_points_and_levels(ParameterInput *pin);



// debug 机制。以后 utils 比较臃肿时可以挪到单独的 debugging.hpp 中定义。


// 通过枚举类型来定义 verbose 的级别，这样做比用数字更加清晰
// 为抑制信息，正值为 debug 信息，0 为正常运行时下应该输出的信息
enum VerboseLevel
{
  //? 目前还没想好 verbose 应该用来输出什么信息
  // 基准是 0 （verbose 的默认值），也就是正常运行时输出：ERROR 和 WARNING
  VERBOSE_NONE = -100,   // 抑制
  VERBOSE_ERROR = -50,   // 遇到严重错误时
  VERBOSE_WARNING = -10, // 可能导致错误的情况，需要检查
  VERBOSE_INFO = 10,     // 一般信息
};
enum DebugLevel
{
  DEBUG_NONE = 0,       // 不输出 debug 信息
  DEBUG_Main = 20,      // 整个主程序的步骤，不包括主循环中的每个时间步
  DEBUG_TimeStep = 30,  // 主循环中每个时间步的信息
  DEBUG_Mesh = 50,      // 调用 Mesh 时相关的信息
  DEBUG_MeshBlock = 60, // 遍历 MeshBlock 相关的信息
  DEBUG_Cell = 80,      // 遍历每个格子相关的信息
  DEBUG_ALL = 100       // 详尽的所有信息
};







#endif // CKJ_CODE_UTILS_HPP_