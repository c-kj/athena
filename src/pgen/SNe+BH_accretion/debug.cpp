// C++ headers
#include <cstdarg>
#include <iostream>
#include <array>
#include <limits>

// Athena++ headers
#include "../../globals.hpp"
#include "../../parameter_input.hpp"
#include "../../eos/eos.hpp"

// 自定义头文件
#include "debug.hpp"
#include "utils.hpp"  // ensure_parent_directory_exists
#include "../ckj_code/ckj_plugin.hpp" // ckj_plugin::current_stage
#include "../SNe+BH_accretion.hpp"  // Abundance


Debug::Debug(Mesh *pmesh, ParameterInput *pin, std::unique_ptr<ProgressReport>& progress_report):
  pmesh(pmesh),
  progress_report(progress_report)
{
  // 从 input file 中读取参数
  level = pin->GetOrAddInteger("debug", "debug", DEBUG_NONE);
  level_on_start = level;
  debug_filepath = pin->GetOrAddString("debug", "debug_filepath", "info/debug_info.txt");
  // verbose_level = pin->GetOrAddInteger("debug", "verbose", 0);
  
  init_debug_stream();  // 不论是否开启 debug，都初始化 debug_stream。否则后续主循环中再开启的话，路径变到 output 下面去了
}


void Debug::init_debug_stream() {
  ensure_parent_directory_exists(debug_filepath);  // 在所有 MPI Rank 上检查/创建目录，不太确定这样是否会在并行时出问题
  if (Globals::my_rank == 0) {
    std::cout << "DEBUG: saving debug info to file " << debug_filepath << std::endl;
  }
  if (not debug_stream.is_open()) { // 如果 debug_stream 没有打开，则打开它
    debug_stream.open(debug_filepath, std::ios::app);
  }
}



void Debug::CheckCell(std::string check_position, MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  Units* punit = pmb->pmy_mesh->punit;
  Real gamma = pmb->peos->GetGamma();
  Real mu = Abundance::mu;

  int ncycle = pmb->pmy_mesh->ncycle;
  static int first_bug_ncycle = -1;
  static int first_negative_ncycle = -1;

  constexpr int var_E_thermal = IPR+1;
  constexpr int var_T = IPR+2;
  constexpr int var_T_div = IPR+3;

  //TEMP 检查 dt。这里用了侵入式写法，临时把 mesh.hpp 中的 private 注释掉了
  std::array<Real, 4> dt_arr {pmb->new_block_dt_, pmb->new_block_dt_hyperbolic_, pmb->new_block_dt_parabolic_, pmb->new_block_dt_user_};
  for (int i=0; i<dt_arr.size(); ++i) {
    Real dt_ = dt_arr[i];
    if (!fpcheck::is_normal(dt_)) {
      debug_stream << "Check dt:"
                  << check_position << ": " 
                  << "ncycle=" << ncycle << ", stage=" << ckj_plugin::current_stage 
                  << ", at gid=" << pmb->gid
                  << ", dt_arr(" << i << ") = " << dt_
                  << std::endl;
    }
  }

  // 遍历每个 cell
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real r = std::sqrt(x*x + y*y + z*z);

        std::stringstream cell_issues;
        bool has_issues = false;
        bool has_non_finite = false;

        for (std::string check: {"non_finite", "negative"}) {  // 遍历检查项
          for (std::string vars_str: {"cons", "prim"}) {

            // 根据 cons 还是 prim 计算 E_thermal、T_cgs 等变量
            auto& vars = vars_str == "cons" ? cons : prim;
            Real rho, P, E_thermal;
            if (vars_str == "prim") {  // 使用 prim 或 cons 来计算 rho、P
              rho = prim(IDN,k,j,i); 
              P = prim(IPR,k,j,i);
              E_thermal = P / (gamma - 1.0);
            } else {
              rho = cons(IDN,k,j,i);
              //TEMP 用 rho 还是 rho_inv ?
              Real rho_inv = 1.0 / rho; 
              E_thermal = cons(IEN,k,j,i) - 0.5 * rho_inv * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));  
              P = (gamma - 1.0) * E_thermal;
            }

            for (int n=0; n <= var_T_div; ++n) { // 遍历每个变量  //TEMP 改为 var_T_div
              Real value;
              if (n <= IPR) {
                value = vars(n,k,j,i);
              } else if (n == var_E_thermal) {
                value = E_thermal;
              } else if (n == var_T) {
                Real T_cgs = P / rho * mu * punit->code_temperature_mu_cgs;
                value = T_cgs;
              } else if (n == var_T_div) { //TEMP 尝试 /rho 和 *rho_inv 有何区别
                Real E_thermal_div = cons(IEN,k,j,i) - 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i))) / rho;
                Real P_div = (gamma - 1.0) * E_thermal_div;
                Real T_cgs_div = P_div / rho * mu * punit->code_temperature_mu_cgs;
                value = T_cgs_div;
                if (vars_str == "prim") { value =  P / rho * mu * punit->code_temperature_mu_cgs; } // 如果是 prim，就还是用 T_cgs 而非 T_cgs_div
              }

              // 触发任何一个条件，则标明这个 cell 有问题，并在 cell_issues 中记录
              if ((check == "non_finite" and !fpcheck::is_finite(value)) or   // non_finite 包括 nan 和 inf
                  (check == "negative" and (n == IDN or n >= IPR) and value <= 0.0)) {  // 其实是检测 non-positive。只对 rho、P、E_thermal、T 等这些量检测
                has_issues = true;
                if (check == "non_finite") { has_non_finite = true; }
                cell_issues << vars_str << "(" << n << ")=" << value 
                            << "[" << check << "]"       // 输出 check 的种类 [non_finite] 还是 [negative]
                            << ", ";
              }
            }
          }
        }

        //TEMP 输出 T_prim_cgs
        Real T_prim_cgs = prim(IPR,k,j,i) / prim(IDN,k,j,i) * mu * punit->code_temperature_mu_cgs;
        if (has_issues and ncycle - first_negative_ncycle < 100 or has_non_finite) { // negative 仅记录前 100 个 cycle，避免太多。每个 MPI rank 分别计数
          debug_stream << check_position << ": " 
                      << "ncycle=" << ncycle << ", stage=" << ckj_plugin::current_stage 
                      << ", at gid=" << pmb->gid
                      << ", (k,j,i)=(" << k << "," << j << "," << i << "), "
                      << cell_issues.str()
                      << "time=" << time 
                      << ", dt=" << dt 
                      << ", (x,y,z,r)=(" << x << "," << y << "," << z << "," << r << "), "
                      << "T_prim_cgs = " << T_prim_cgs
                      << std::endl;

        }
        
        if (has_non_finite) { // 只有 cell 中含有 nan_finite 的情况下才认为是 bug，记录 ncycle。如果一直只有 negative 但没有 nan，目前认为不算 bug。
          if (first_bug_ncycle == -1) { first_bug_ncycle = ncycle; }
          if (ncycle - first_bug_ncycle > 3) {
            throw std::runtime_error("在连续 3 个 cycle 中，检测到 cell 中有 NaN，终止模拟！\n");
          }
        }
        if (has_issues and first_negative_ncycle == -1) { first_negative_ncycle = ncycle; }
        
      }
    }
  }
}



void Debug::check_zero_dt() {
  static int zero_dt_count = 0;  // 连续出现 dt == 0 的 cycle 数
  static int level_original = level;  // 被本函数修改前的 level

  if (pmesh->dt == 0.0) { // 如果出现 dt == 0 的 bug
    if (zero_dt_count == 0) {  // 初次检测到
      if (Globals::my_rank == 0) {
        std::cerr << "在 ncycle = " << pmesh->ncycle << " 时开始检测到 dt==0，开启 debug = DEBUG_ALL。\n";
      }
      level = DEBUG_ALL;  // 开启 debug flag（在每个 MPI rank 上）
    }
    zero_dt_count++;
    if (zero_dt_count > 100 and Globals::my_rank == 0) {
      std::string msg = "在连续 100 个 cycle 中，dt 都为 0，终止模拟！";
      if (progress_report) { progress_report->stream << msg << std::endl; } // 先判空，防止 seg fault
      throw std::runtime_error(msg);
    }
  } else { // 正常情况
    if (zero_dt_count > 0) {   // 如果 dt 从 0 恢复正常了
      zero_dt_count = 0;       // 清零计数
      level = level_original;  // 恢复原来的 debug level
    }
  }
}



void Debug::close() {
  if (debug_stream.is_open()) {
    debug_stream << "debug_stream is open, closing \n" << std::endl;
    debug_stream.close();
  }
}


// 为全局的 debug 指针提供默认定义，从而在其他文件中可以直接使用（当然，要先初始化）
std::unique_ptr<Debug> debug;




#if false // deprecated

// 自定义的用于将 debug 信息同时 print 和输出到文件的函数
// 定义一个函数，它接收一个输出流、一个字符串前缀、一个格式化字符串和一些可变参数
//! 注意！不知道这个函数如果用于 MPI 并行计算时，会不会出问题！可能会抢占式写入？？？最好是只用在本地 debug 时
//! 对于 MeshBlock 级别以上，OpenMP or MPI 时也是危险的，同时往一个流里面写入
void printf_and_save_to_stream(std::ostream &stream, const char *format, ...) {
  va_list args;           // 定义一个 va_list 类型的变量，用于存储可变参数列表
  va_start(args, format); // 使用 va_start 函数初始化 args，使其包含从 format 之后的所有参数

  char buffer[256]; // 定义一个字符数组，用于存储格式化后的字符串

  // 使用 vsnprintf 函数将格式化的字符串写入 buffer，vsnprintf 函数会根据 format 和 args 生成一个字符串，并确保不会超过 buffer 的大小
  vsnprintf(buffer, sizeof(buffer), format, args);

  printf("%s", buffer); // 使用 printf 函数将 buffer 输出到屏幕

  // 检查 stream 是否有效
  if (stream)
  {
      // 如果 stream 有效，将 output 输出到 stream
      stream << buffer;
  }
  else
  {
      // 如果 stream 无效，输出一条错误消息
      printf("stream is not open !!!!! \n");
  }

  // 使用 va_end 函数结束可变参数列表
  va_end(args); 
}

#endif