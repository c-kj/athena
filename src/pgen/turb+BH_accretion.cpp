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
  // 自定义 AMR 机制所使用的变量
  Real time_start_AMR;
  int my_root_level;     // Mesh 的 root_level 是 private 的，所以需要自己定义一个变量来存储它
  std::vector<std::array<Real, 3>> point_list;
  std::vector<int> level_list;

  int verbose, debug;
  std::string debug_filepath, verbose_filepath;
  std::ofstream debug_stream;

  std::string init_cond_type;
  Real rho_at_boundary, E_tot_at_boundary, power_law_index, alpha;
  Real approx_Bondi_rho_profile(Real alpha, Real R_Bondi, Real r) {
    return pow(1 + alpha * R_Bondi / r, 1.5);
  }

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
  //! 对于 MeshBlock 级别以上，OpenMP or MPI 时也是危险的，同时往一个流里面写入
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

  // 用于从 input file 中读取自定义的 AMR 参数，以 point_1 = 0, 1, 2.0 这种形式给出点的坐标
  std::tuple<std::vector<std::array<Real, 3>>, std::vector<int>> get_AMR_points_and_levels(ParameterInput *pin) {
    std::string block_name, point_name, level_name, num;
    std::array<Real, 3> point;
    std::vector<std::array<Real, 3>> point_list;
    std::vector<int> level_list;
    int level, default_level;

    default_level = pin->GetInteger("mesh", "numlevel") - 1; // numlevel 是从 1 开始的，所以要减 1
    block_name = "AMR/point";                    // 从 input file 中读取的 Block 的名字

    for (int i = 1; true; i++) {                // 让 i 递增，直到找不到 point_i 为止
      point_name = "point_" + std::to_string(i);
      level_name = "level_" + std::to_string(i);
      if (pin->DoesParameterExist(block_name, point_name)) {
        std::stringstream ss(pin->GetString(block_name, point_name));
        // 把逗号分隔的三个数字读入 point 数组
        for (int j = 0; j < 3; j++) {
          std::getline(ss, num, ',');
          point[j] = std::stod(num);
        }
        level = pin->GetOrAddInteger(block_name, level_name, default_level);   // 把 level_i 读入 level
        // 把 point 和 level 添加到 point_list 和 level_list 中
        point_list.push_back(point);
        level_list.push_back(level);

      } else {
        break;
      }
    }
    return std::make_tuple(point_list, level_list);    // 返回一个 tuple
  }

}
using namespace my_utils;


// 通过判断 MeshBlock 是否包含指定的点，确定是否需要细化
int RefinementCondition_Point(MeshBlock *pmb, std::array<Real, 3> point, int level_limit=100000) {
  RegionSize size = pmb->block_size;
  int current_level = pmb->loc.level - my_root_level;
  if (current_level >= level_limit) return 0; // 如果当前 level 大于等于 level_limit，则不再细化

  Real x1 = point[0], x2 = point[1], x3 = point[2];
  if (size.x1min <= x1 && x1 <= size.x1max &&
      size.x2min <= x2 && x2 <= size.x2max &&
      size.x3min <= x3 && x3 <= size.x3max) {
    return 1;
  }
  return 0;
}


//? 是否有必要把函数的声明和定义分开？声明放在前面，定义放在后面？
// 测试能不能用 AMR 来充当 SMR，初始化湍流
int RefinementCondition(MeshBlock *pmb) {
  // 在 time_start_AMR 之前，都不进行 AMR，直接返回 0
  Real time = pmb->pmy_mesh->time;
  if (time < time_start_AMR) return 0;

  // debug 消息，放在前面，免得还没输出就 return 了
  if (debug >= DEBUG_Mesh && pmb->gid == 0) {
    printf_and_save_to_stream(debug_stream, "DEBUG_Mesh: RefinementCondition is called at time = %f \n", time);
  }

  //* 这里涉及有些复杂的逻辑判断：每一个 condition 返回的是 -1,0,1 之一，如果有 -1 和 1 的冲突那么需要仔细考量优先级。
  for (int i = 0; i < point_list.size(); i++) {
    if (RefinementCondition_Point(pmb, point_list[i], level_list[i]) == 1) return 1; // 如果 MeshBlock 包含任何一个指定的点，则进行细化
  }
  // if (RefinementCondition_Point(pmb, {0,0,0}) != 0) return 1; // 如果包含原点，则进行细化
  // 这里这么写是为了后续扩展
  
  return 0;
}



void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

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
  
  // 根据 debug 信息，自定义源项在每个时间步会被 call 两次：一次是在时间步开始时，另一次在时间步的一半处。可能这就是 operator splitting？
  if (debug >= DEBUG_Mesh && pmb->gid == 0 ){
    printf_and_save_to_stream(debug_stream, "DEBUG_Mesh: Source Term is called at time = %f \n", time);
  }

  return;
}





//========================================================================================
// 以下是 Athena++ 提供的接口
//========================================================================================

// 调用时机：程序启动时（包括 restart 时）
// 函数用途：初始化 Mesh 上的全局变量（各个 MeshBlock 共享）；enroll 各种自定义函数
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

  my_root_level = root_level;     // Mesh 的 root_level 是 private 的，所以需要自己定义一个变量来存储它

  // 从 input file 中读取参数，存到 pgen 文件的全局变量中
  GM_BH = pin->GetReal("problem","GM_BH");
  R_in = pin->GetReal("problem","R_in");
  R_out = pin->GetReal("problem","R_out");
  rho_in_BH = pin->GetReal("problem","rho_in_BH");
  rho_init = pin->GetReal("problem","rho_init");
  E_tot_init = pin->GetReal("problem", "E_tot_init");

  // 读取自定义的 AMR 参数
  if (adaptive) {
    time_start_AMR = pin->GetOrAddReal("AMR", "time_start_AMR", 0.0);

    std::tie(point_list, level_list) = get_AMR_points_and_levels(pin);  // 用 tie 函数将返回的 tuple 解包
  }

  // <debug> 里的参数
  debug = pin->GetOrAddInteger("debug", "debug", 0);
  if (debug > DEBUG_NONE) {
    debug_filepath = pin->GetOrAddString("debug", "debug_filepath", "info/debug_info.txt");
    ensure_parent_directory_exists(debug_filepath);
    debug_stream.open(debug_filepath, std::ios::app); // 用于输出 debug 信息的文件流
    printf("DEBUG: saving debug info to file %s \n", debug_filepath.c_str());
  }

  // 将自定义的源项注册到 Athena++ 中
  // 调用时机：在每个时间步中调用两次，一次在开头，一次在中间
  EnrollUserExplicitSourceFunction(SMBH_grav);

  // 将自定义的 AMR 条件注册到 Athena++ 中
  if(adaptive){
    EnrollUserRefinementCondition(RefinementCondition);
  }

  return;
}

// 调用时机：程序启动时（包括 restart 时）
// 函数用途：初始化每个 MeshBlock 上的局部变量（包括数组数据）
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // 之前我用来传递 input 参数，后来改成了用 pgen 的全局变量
  return;
}


// 调用时机：程序启动时（应该不包括 restart 时？）
// 函数用途：设定初始条件（每个 MeshBlock 的）
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  init_cond_type = pin->GetOrAddString("initial_condition","init_cond_type","uniform");
  power_law_index = pin->GetOrAddReal("initial_condition","power_law_index",0.0);   // 一般情况下是负数
  alpha = pin->GetOrAddReal("initial_condition","alpha",0.3333333);

  if (init_cond_type == "uniform") {
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

  // 幂律初值，密度和能量遵循相同的幂律，初速度为 0。初始声速为 gamma*(gamma-1)*E_tot/rho
  if (init_cond_type == "power_law") {
    Real x1,x2,x3,r;
    for (int k=ks; k<=ke; k++) {
      x3 = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
        x2 = pcoord->x2v(j);
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          r = sqrt(x1*x1 + x2*x2 + x3*x3);
          phydro->u(IDN,k,j,i) = rho_init * pow(r/R_out, power_law_index);

          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = E_tot_init * pow(r/R_out, power_law_index);
          }
        }
      }
    }
  }

  // approximate_Bondi: 初速度 = 0，密度和能量遵循 Bondi profile 的近似解
  if (init_cond_type == "approximate_Bondi") {
    Real gamma = peos->GetGamma();
    Real R_Bondi = 2 * GM_BH / (gamma*(gamma-1)*E_tot_init/rho_init);
    //TODO 目前这里是把 init 的值直接理解为边界值，然后换算出无穷远的值。以后可以考虑修改
    Real rho_at_boundary = rho_init;
    Real E_tot_at_boundary = E_tot_init;  
    Real rho_inf = rho_at_boundary / approx_Bondi_rho_profile(alpha, R_Bondi, R_out);


    Real x1,x2,x3,r;
    for (int k=ks; k<=ke; k++) {
      x3 = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
        x2 = pcoord->x2v(j);
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          r = sqrt(x1*x1 + x2*x2 + x3*x3);
          phydro->u(IDN,k,j,i) = rho_inf * approx_Bondi_rho_profile(alpha, R_Bondi, r);

          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = E_tot_at_boundary * pow(phydro->u(IDN,k,j,i) / rho_at_boundary, gamma) ;
          }
        }
      }
    }
  }

}


// 调用时机：在每个时间步的结尾。每个 MeshBlock 调用一次
// 函数用途：仅用于分析，不用于操作数据
void MeshBlock::UserWorkInLoop() {
  return;
}


// 调用时机：在每个时间步的结尾。每个节点调用一次
// 函数用途：可用于从 MeshBlocks 收集结果；使用 MPI all-to-all 通信
void Mesh::UserWorkInLoop() {
  // 在时间步内调用 Mesh::UserWorkInLoop 时，步长内主要的演化已经完成，但 time 还未更新。
  // 所以「时间步结尾的时间」要加上 dt
  Real end_time = time + dt; 
  if (debug >= DEBUG_TimeStep) {
    printf_and_save_to_stream(debug_stream, "DEBUG_TimeStep: Mesh::UserWorkInLoop is called at end_time = %f \n", end_time);
  }
  return;
}


// 调用时机：每个将要输出 output 的时间步的末尾。和 Mesh::UserWorkInLoop 谁先？
// 函数用途：计算用户定义的输出变量
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  return;
}


// 调用时机：模拟结束时
// 函数用途：清理 MeshBlock::InitUserMeshData 分配的资源
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
