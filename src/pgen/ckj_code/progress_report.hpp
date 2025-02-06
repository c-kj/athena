#pragma once
// 用于进行 progress report 的类。输出模拟进度、耗时等信息

#include <string>   // std::string
#include <fstream>  // std::ofstream
#include <chrono>   // std::chrono
#include <sstream>  // std::stringstream
#include <iomanip>  // std::fixed, std::setprecision

// Athena++ headers
#include "../../athena.hpp"
#include "../../globals.hpp"
#include "../../mesh/mesh.hpp"

// 自定义的头文件（各种模块）
#include "utils.hpp"  // ensure_parent_directory_exists, format_duration


// 用于进行 progress report 的类。输出模拟进度、耗时等信息
// 在内部判断 MPI rank，只在 rank 0 上打开文件流、执行输出
struct ProgressReport {
  // 构造函数。在 Mesh::InitUserMeshData 中构造一个实例
  ProgressReport(Mesh *pmesh, ParameterInput *pin);

  // 判断当前 cycle 内是否要 report
  bool is_time_to_report() const; 
  // 执行汇报操作
  void report(bool final=false);  // final 为 true 时，表示最后一次 report，会输出总结信息
  
  //* 以下 time 指模拟中的 code time 的时间 pmesh->time，而 t 指现实中的 wall time
  const Mesh* const pmesh;
  std::ofstream stream;
  const int dcycle_report;  // 每隔多少个 cycle 输出一次 hst
  const std::chrono::high_resolution_clock::time_point t_start;  // 初始化时的时间点

  struct { // 上次 report 时的信息
    std::chrono::high_resolution_clock::time_point t;
    Real time;   // 用于计算这个汇报周期的耗时
    int ncycle;  // 用于确定下次汇报是什么时候
  } last_report; 
};


// 以下是函数的实现。因为不多，就不单独写到 .cpp 文件中了
// 在头文件中定义时，需要用 inline 避免多重定义错误


inline ProgressReport::ProgressReport(Mesh *pmesh, ParameterInput *pin):
    pmesh(pmesh),
    dcycle_report(pin->GetOrAddInteger("debug","dcycle_report", 1000)),
    t_start(std::chrono::high_resolution_clock::now()) 
{
  if (Globals::my_rank == 0) {
    // 初始化 last_report。这里不直接使用初始化列表，是为了让名称对应得更加清晰
    last_report.t = t_start;
    last_report.time = pmesh->time;
    last_report.ncycle = pmesh->ncycle;
  
    // progress_report 的输出文件流。只在 rank 0 上打开（注意：在其他 rank 上并不打开）。
    std::string filepath = "info/progress_report.txt";
    ensure_parent_directory_exists(filepath);
    stream.open(filepath, std::ios::app); // 用于输出 progress report 的文件流。使用 append 模式是考虑到有可能是 restart，这时不覆盖之前的内容。
    if (!stream.is_open()) {  // 检查是否成功打开文件
      throw std::runtime_error("Failed to open file: " + filepath);
    }
  }
}


inline bool ProgressReport::is_time_to_report() const {
  int ncycle = pmesh->ncycle;
  int dcycle = (ncycle <= dcycle_report) ? (dcycle_report/10) : dcycle_report;  // 实际使用的 dcycle_report，在刚开始的一段时间内输出更频繁一些。这里用 int 的除法自动向下取整

  //FUTURE 可以不仅用 dcycle 控制，还有现实中的时长间隔
  return ( ncycle - last_report.ncycle >= dcycle );  // 不采用 % 取余的方法，因为有可能 Athena++ 的 ncycle_out 不为 1，不是每个 cycle 都会有输出 
}


inline void ProgressReport::report(bool final) {
  // 进度汇报、耗时估算
  // 预报：预期总时长、预期剩余时长。多种估算方法：历史、瞬时

  if (Globals::my_rank == 0) {
    int ncycle = pmesh->ncycle;
    Real time  = pmesh->time;
    Real dt    = pmesh->dt;
    Real tlim  = pmesh->tlim;
    Real nlim  = pmesh->nlim;
    Real start_time = pmesh->start_time;
    Units *punit = pmesh->punit;

    using namespace std::chrono;

    auto t_now = high_resolution_clock::now();
    duration<double> elapsed = t_now - t_start;
    duration<double> elapsed_since_last_report = t_now - last_report.t;

    Real total_time = (tlim - start_time);  // 模拟总时长，单位为 code time
    

    // 进度汇报
    std::stringstream text;  // 当前 cycle 输出的 progress report 文本。用于同时输出到文件和 std::cout
    text << "---------------------------------------------------------------\n";
    text << std::fixed << std::setprecision(2); // 百分比使用固定小数点
    text << "cycle = " << ncycle << ". 当前耗时: " <<  format_duration(elapsed) << ", 进度 " << time/total_time * 100 << "% \n";
    text << std::scientific << std::setprecision(2);  // 以下数字使用科学计数法。时分秒除外
    text << "当前模拟中 time = " << time << " (" << time * punit->million_yr_code << " Myr), 全长 " << total_time << " (" << total_time * punit->million_yr_code << " Myr)\n";
    text << "最近 " << ncycle - last_report.ncycle << " 个 cycle 耗时: " << format_duration(elapsed_since_last_report) << ", 模拟中时间推进了 " << (time - last_report.time) << " code_time (" << (time - last_report.time) * punit->million_yr_code << " Myr)\n";

    // 耗时估算
    duration<double> t_total_all = elapsed / time * total_time;  // 基于历史耗时，估计总时长
    duration<double> t_remains_all = t_total_all - elapsed;               // 根据总时长估计剩余时长
    text << "根据历史耗时，预计总时长: " << format_duration(t_total_all)     << "，预计剩余时长: " << format_duration(t_remains_all) << "\n";
    duration<double> t_remains_recent = elapsed_since_last_report / (time - last_report.time) * (tlim - time);   // 根据最近耗时，估计剩余时长
    duration<double> t_total_recent   = elapsed + t_remains_recent;                                        // 根据剩余时长，估计总时长
    text << "根据最近耗时，预计总时长: " << format_duration(t_total_recent) << "，预计剩余时长: " << format_duration(t_remains_recent) << "\n";
    text << "---------------------------------------------------------------\n";

    if (final) {
      text << "模拟结束。\n";
    }

    // 把 progress report 同时输出到单独的文件 progress_report_stream 和 std::cout 中
    stream << text.str();
    stream.flush();  // 刷新输出流，这样才会立刻写入文件
    std::cout << text.str();

    // 更新 last_report
    last_report.t = t_now;
    last_report.time = time;
    last_report.ncycle = ncycle;
  }
}
