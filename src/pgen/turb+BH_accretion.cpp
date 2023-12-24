//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//! \brief Problem generator for turbulence driver

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//========================================================================================
// 自定义的变量、函数等
//========================================================================================

// 自定义的头文件导入
#include <fstream>  // std::ofstream
#include <iostream> // std::cout, std::endl
#include <cstdarg>  // va_list, va_start, va_end, vsprintf 等用于处理变长参数列表
#include <sys/stat.h>  // 为了使用 stat 和 mkdir 函数
// #include <string>  // 这里无需 include，因为在 parameter_input.hpp 中已经 include 了


// 声明自定义的全局变量，从 input file 中读取
// 这些变量的赋值是在 Mesh::InitUserMeshData 中完成的
// 使用 namespace 创建一个命名空间 my_vars，然后再 using namespace 导入它
// 这样做的好处是，所有自定义的变量都放在一起，在大纲层级中更加清晰
namespace my_vars {
  Real GM_BH, R_in, R_out, rho_in_BH, rho_init, E_tot_init;
  Real current_time, time_start_AMR;
  int verbose, debug;
  std::string debug_filepath, verbose_filepath;
  std::ofstream debug_stream;

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
}
using namespace my_vars;

// 自定义的工具函数
namespace my_utils {
  // 自定义的用于将 debug 信息同时 print 和输出到文件的函数
  // 定义一个函数，它接收一个输出流、一个字符串前缀、一个格式化字符串和一些可变参数
  //! 注意！不知道这个函数如果用于 MPI 并行计算时，会不会出问题！可能会抢占式写入？？？最好是只用在本地 debug 时
  void printf_and_save_to_stream(std::ostream& stream, const char* format, ...) {
    va_list args;  // 定义一个 va_list 类型的变量，用于存储可变参数列表
    va_start(args, format); // 使用 va_start 函数初始化 args，使其包含从 format 之后的所有参数

    char buffer[256];  // 定义一个字符数组，用于存储格式化后的字符串

    // 使用 vsnprintf 函数将格式化的字符串写入 buffer，vsnprintf 函数会根据 format 和 args 生成一个字符串，并确保不会超过 buffer 的大小
    vsnprintf(buffer, sizeof(buffer), format, args); 

    printf("%s", buffer);   // 使用 printf 函数将 buffer 输出到屏幕

    // 检查 stream 是否有效
    if (stream) {
      // 如果 stream 有效，将 output 输出到 stream
      stream << buffer; 
    } else {
      // 如果 stream 无效，输出一条错误消息
      printf("stream is not open !!!!! \n");
    }

    // 使用 va_end 函数结束可变参数列表
    va_end(args); 
  }


  // 检查并创建目录的函数
  void check_and_create_directory(const std::string& path) {
    struct stat info;  // 存储文件或目录的信息
    if (stat(path.c_str(), &info) != 0) {
      mkdir(path.c_str(), 0755);  // 如果路径不存在，创建它
    } else if (!(info.st_mode & S_IFDIR)) {
      throw std::runtime_error(path + " exists but is not a directory");  // 如果路径存在，但不是一个目录，抛出错误
    }
  }

  // 检查一个路径的父目录是否存在，如果不存在则创建它。
  // 会创建 path 中最后一个 / 之前的全部路径，因此结尾是否有 / 是不同的
  void ensure_parent_directory_exists(const std::string& path) {
    std::string partial_path;  // 存储路径的部分内容

    for (char c : path) {
      partial_path += c;  // 将当前字符添加到部分路径中
      if (c == '/') {
        check_and_create_directory(partial_path);  // 检查并创建部分路径
      }
    }
  }

}
using namespace my_utils;


//? 是否有必要把函数的声明和定义分开？声明放在前面，定义放在后面？
// 测试能不能用 AMR 来充当 SMR，初始化湍流
int RefinementCondition(MeshBlock *pmb) {
  // 在 time_start_AMR 之前，都不进行 AMR，直接返回 0
  if (current_time < time_start_AMR) return 0;

  Real x,y,z,dx,dy,dz;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    dz = pmb->pcoord->dx3f(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      dy = pmb->pcoord->dx2f(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        dx = pmb->pcoord->dx1f(i);

        // 如果当前 MeshBlock 中有任何一个格点，其 xyz 坐标绝对值均小于其格子边长，则认定该 MeshBlock 紧贴坐标原点，需要细化
        if (std::abs(x) < dx && std::abs(y) < dy && std::abs(z) < dz) return 1;
      }
    }
  }
  
  if (debug >= DEBUG_Mesh && pmb->gid == 0) {
    printf_and_save_to_stream(debug_stream, "DEBUG_Mesh: RefinementCondition called at time = %f \n", current_time);
  }

  return 0;
}



void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  current_time = time;

  Real x,y,z,r,r3,vx,vy,vz; 
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        r = sqrt(x*x + y*y + z*z);
        r3 = r*r*r;
        
        vx = prim(IVX,k,j,i);
        vy = prim(IVY,k,j,i);
        vz = prim(IVZ,k,j,i);
        Real& rho = cons(IDN,k,j,i); // rho 是引用，修改 rho 会修改 cons(IDN,k,j,i) 的值

        // TODO 只对黑洞外的区域施加引力，这样处理正确吗？
        // TODO 这里需要 && r <= R_out 吗？
        if (r >= R_in ) {
          cons(IM1,k,j,i) += - GM_BH*x/r3 * rho * dt;
          cons(IM2,k,j,i) += - GM_BH*y/r3 * rho * dt;
          cons(IM3,k,j,i) += - GM_BH*z/r3 * rho * dt;

          // 能量的改变
          // ? 我看 Athena++ 的源代码都是判断了 if NON_BAROTROPIC_EOS，但这玩意默认是 True，不知道啥意思… 绝热方程是不是 barotropic？从模拟结果上来看好像不是
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) += - GM_BH*(x*vx + y*vy + z*vz)/r3 * rho * dt; // 这里可以稍微优化，快一点
          }

        }

        // 进入黑洞的处理
        if (r < R_in) {
          // rho = std::min(rho_floor, rho); // 进入黑洞后，密度取 min(rho_floor, rho) // TODO 这个正确吗？有必要吗？可能自带约束？
          rho = rho_in_BH;

          // TODO 修改动量使其不外流。这样做正确吗？
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;

          // 修改能量
          cons(IEN,k,j,i) = 0.0 + 0.0; // TODO 应该怎么修改能量？尤其是内能怎么计算？ shocks?

        }

        // TODO 外部区域的处理，应该怎么做？
        // 目前是设定为保持初始值
        // TODO 或许可以改为继承该格子之前的值？但怎么操作呢？在 source term 里面不知道能不能做这件事
        if (r > R_out) {
          rho = rho_init;
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          cons(IEN,k,j,i) = E_tot_init;
        }

        // debug 信息：逐个格点检查 rho 是否为 nan
        if (debug >= DEBUG_Cell && std::isnan(cons(IDN,k,j,i)) ) {
          printf_and_save_to_stream(debug_stream, "DEBUG_Cell: rho is nan at x,y,z,r = (%f, %f, %f, %f) \n", x, y, z, r);
        }

      }
    }
  }
  
  if (debug >= DEBUG_Mesh && pmb->gid == 0 ){
    printf("saving debug info to file %s \n", debug_filepath.c_str());
    printf_and_save_to_stream(debug_stream, "DEBUG_Mesh: calling Source Term at time = %f \n", time);
  }

  return;
}





//========================================================================================
// 以下是 Athena++ 提供的接口
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem","four_pi_G");
    SetFourPiG(four_pi_G);
  }
  
  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for impulsively driven turbulence
  // turb_flag = 3 for continuously driven turbulence
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  // 从 input file 中读取参数，存到 pgen 文件的全局变量中
  GM_BH = pin->GetReal("problem","GM_BH");
  R_in = pin->GetReal("problem","R_in");
  R_out = pin->GetReal("problem","R_out");
  rho_in_BH = pin->GetReal("problem","rho_in_BH");
  rho_init = pin->GetReal("problem","rho_init");
  E_tot_init = pin->GetReal("problem", "E_tot_init");

  time_start_AMR = pin->GetOrAddReal("problem", "time_start_AMR", 0.0);

  // <debug> 里的参数
  // verbose = pin->GetOrAddInteger("debug", "verbose", 0);
  // verbose_filepath = pin->GetOrAddString("debug", "verbose_filepath", "info/verbose_info.txt");
  debug = pin->GetOrAddInteger("debug", "debug", 0);
  if (debug >= DEBUG_NONE) {
    debug_filepath = pin->GetOrAddString("debug", "debug_filepath", "info/debug_info.txt");
    ensure_parent_directory_exists(debug_filepath);
    debug_stream.open(debug_filepath, std::ios::app); // 用于输出 debug 信息的文件流
  }

  // 将自定义的源项注册到 Athena++ 中
  EnrollUserExplicitSourceFunction(SMBH_grav);

  // 将自定义的 AMR 条件注册到 Athena++ 中
  if(adaptive){
    EnrollUserRefinementCondition(RefinementCondition);
  }

  return;
}


// 这里可以设定属于每个 MeshBlock 的自定义数组数据。之前我用来传递 input 参数，后来改成了用 pgen 的全局变量
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  return;
}



// MeshBlock::ProblemGenerator 用于设定初始条件（每个 MeshBlock 的）
// defined in either the prob file or default_pgen.cpp in ../pgen/
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho_init;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = E_tot_init;
        }
      }
    }
  }
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // 这个函数用于在模拟结束后做一些事情
  // TODO
  if (debug >= DEBUG_Main) {
    printf_and_save_to_stream(debug_stream, "DEBUG_Main: Mesh::UserWorkAfterLoop is called \n");
  }

  // 关闭 debug_stream
  if (debug_stream.is_open()) { 
    printf("debug_stream is open, closing \n");
    debug_stream.close();
  }
}
