#pragma once


// C++ headers
#include <memory> // std::unique_ptr

//TODO 检查这些 include 有无必要?
// #include <fstream>  
// #include <string>
// #include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"

// 自定义的头文件
#include "progress_report.hpp" // progress_report 指针


class Debug {
public:
  Debug(Mesh *pmesh, ParameterInput *pin, std::unique_ptr<ProgressReport>& progress_report);

  // 初始化 debug_stream：确保文件存在、确保 stream 打开
  void init_debug_stream();

  // 检查单元格的值是否有问题
  void CheckCell(std::string check_position, MeshBlock *pmb, const Real time, const Real dt,
                 const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                 const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                 AthenaArray<Real> &cons_scalar);
                 
  // 检查 dt == 0 的 bug，及时退出。目前出现在本该 Crash 但加了 floor 导致一直停滞，模拟无法结束的情形
  void check_zero_dt();

  //FUTURE 即便 debug level 为 0 的情况下，依然粗略地检查是否有 nan，避免浪费金钱。但怎么做才能让代价较小？

  // 关闭 debug_stream
  void close();


  Mesh *pmesh;
  std::unique_ptr<ProgressReport>& progress_report;  // 需要保留 progress_report 智能指针的引用，从而向 progress_report 的 stream 中输出消息
  // 这里好像有点蠢，不应该用 unique_ptr 的引用？但这样还挺省事的，免得用指针还要担心空指针。暂时也没什么 bug。

  int level;   // debug 级别。在每个 MPI rank 上独立、可变
  int level_on_start;  // 初始时 level 的值。不过 restart 时不会保持以前的值，而是会重新从 pin 读取
  int verbose_level;     // 详细信息级别（目前没有启用）
  std::string debug_filepath;
  std::ofstream debug_stream;  // 用于输出 debug 信息的文件流
};


// 全局的 Debug 机制指针。在 debug.cpp 中定义。其他文件可以直接使用这个指针。
// 需要在 Mesh::InitUserMeshData 中赋值初始化。
// 这个类的设计上应该是单例（在每个 MPI rank 上），但并不做任何检查来保证。
extern std::unique_ptr<Debug> debug;

// deprecated: 这个函数目前不再有调用（全都注释掉了），以后也最好别调用了。出于参考价值放在这里。
void printf_and_save_to_stream(std::ostream &stream, const char *format, ...);




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

// Debug 的级别，通过 debug->level >= xxx 来判断是否需要启用某些 debug 代码段。
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