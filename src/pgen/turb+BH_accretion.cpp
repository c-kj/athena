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

// 自定义函数所需的标准库导入
// #include <fstream>  // std::ofstream
// #include <string>  // 这里无需 include，因为在 parameter_input.hpp 中已经 include 了

// 自定义的头文件
#include "turb+BH_accretion.hpp"
#include "ckj_code/utils.hpp"
#include "ckj_code/region.hpp"
#include "ckj_code/refinement_condition.hpp"
#include "ckj_code/supernova.hpp"
#include "ckj_code/cooling.hpp"


Real approx_Bondi_rho_profile(Real alpha, Real R_Bondi, Real r) {
  return pow(1 + alpha * R_Bondi / r, 1.5);
}

// Source Term 源项
//* 应该只依赖于 prim 而修改 cons。不能修改 prim。不应该依赖于 cons（尤其是因为 cons 可能被前面的源项修改）。
//? 以后可能考虑把各个源项拆分。这样需要反复进入循环好几次，但是更加清晰。也可以考虑把最内层循环内的各个步骤抽象成函数？
void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  //* CoolingSourceTerm 目前放在 SourceTerm 中，以后可能考虑改到 UserWorkInLoop 中
  if (cooling.cooling_flag) {
    cooling.CoolingSourceTerm(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  const Real& mesh_time = pmb->pmy_mesh->time;
  const Real& mesh_dt = pmb->pmy_mesh->dt;
  const Real dt_ratio = dt/mesh_dt;

  auto punit = pmb->pmy_mesh->punit;
  Real SN_energy_unit, SN_mass_unit;
  if (true) { // 这里暂时没有写开关：什么情况下把 SN 的 input 单位解读为 1e51 erg 和 Msun。若有需求再改
    SN_energy_unit = punit->bethe_code;
    SN_mass_unit = punit->solar_mass_code;
  };

  Real inject_ratio = 0.0; // 用于控制超新星能量、质量、动量等的注入比例。对于 vl2 和 rk2 两种方法不是很有必要，因为最优都是只在某一步注入 1.0。

  //* 处理适应于 integrator 的超新星注入时机
  bool SN_flag_in_this_step = false;
  if (SN_flag > 0) {
    // 判断是否应该在这次调用 source term 时注入超新星
    if (integrator == "vl2") { // VL2: 只在校正步注入
      if (dt_ratio == 1.) {
        inject_ratio = 1.;
        SN_flag_in_this_step = true;
      } else if (dt_ratio == .5) {
        inject_ratio = 0.; // 奇怪，这里设为负数都行，不会影响结果。但 >0.25 左右就会炸步长
      } 
    } else if (integrator == "rk2") { // RK2: 只在完整步长的那一步注入 
      if (dt_ratio == 1.) {
        inject_ratio = 0.; // 这里似乎会影响注入的能量/质量，但数值只有一半的贡献。具体为什么暂时不管了。会炸步长。
      } else if (dt_ratio == .5) {
        inject_ratio = 1.;
        SN_flag_in_this_step = true;
      } 
    } else { //* 其他 integrator，目前不处理，让他炸步长吧
      // throw std::runtime_error("Unsupported integrator! Only VL2 and RK2 are supported.");
      inject_ratio = 1.;
      SN_flag_in_this_step = true;
    }
  }

  if ( SN_flag_in_this_step ) { // 只有这一步有可能注入超新星时，才找出需要注入的超新星
    // 找到当前时间步内需要注入的超新星，放入 supernova_to_inject 列表
    supernova_to_inject = {}; // 每次都要先清空
    for (SuperNova& SN : supernova_list) {
      for (Real& SN_time : SN.time_list) {
        if (mesh_time <= SN_time && SN_time < mesh_time + mesh_dt) {
          supernova_to_inject.push_back(&SN); // &SN 解引用，取 SN 的地址，也就是一个指向 SN 的指针
        }
      }
    }
  }

  Real x,y,z,r,r3,vx,vy,vz,rho; 
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
#pragma omp simd  // OpenMP 的 SIMD 并行。在纯 MPI 并行时可能不起作用？
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        r = sqrt(x*x + y*y + z*z);
        r3 = r*r*r;
        
        rho = prim(IDN,k,j,i); 
        vx = prim(IVX,k,j,i);
        vy = prim(IVY,k,j,i);
        vz = prim(IVZ,k,j,i);

        //? 只对黑洞外的区域施加引力，这样处理正确吗？不过在这里对 R_in 以内施加引力没用，接下来会被重设
        //? 这里需要 && r <= R_out 吗？
        if (r >= R_in ) {
          cons(IM1,k,j,i) += - GM_BH*x/r3 * rho * dt;
          cons(IM2,k,j,i) += - GM_BH*y/r3 * rho * dt;
          cons(IM3,k,j,i) += - GM_BH*z/r3 * rho * dt;

          // 能量的改变：非全局正压时需要更新能量密度
          if (NON_BAROTROPIC_EOS) {
            // 这里可以稍微优化，快一点
            cons(IEN,k,j,i) += - GM_BH*(x*vx + y*vy + z*vz)/r3 * rho * dt; // 根据动量的改变，相应地改变动能。内能目前不变。
            //? 这里是否需要考虑把平方项补偿上去？有待测试。
          }

          // 超新星爆炸的能量和质量注入。//* 目前没有实现动量注入，也没有考虑 SN 的运动速度
          if ( SN_flag_in_this_step ) {
            for (SuperNova* SN : supernova_to_inject) {
              if ( SN->energy_region->contains({x,y,z}) ) {
                cons(IEN,k,j,i) += SN->energy_density * inject_ratio * SN_energy_unit;
              }
              if ( SN->mass_region->contains({x,y,z}) ) {
                cons(IDN,k,j,i) += SN->mass_density * inject_ratio * SN_mass_unit;
              }
              // 这里的 debug 信息不再那么有用，将来可以考虑删去
              // if (debug >= DEBUG_Main && pmb->gid == 0 && SN_flag > 0){
              //   if(mesh_time <= SN_time && SN_time < mesh_time + mesh_dt && dt == mesh_dt ) {
                  // printf_and_save_to_stream(debug_stream, 
                  // "DEBUG_Mesh: SN explosion at time = %f, mesh_time = %f, dt = %f, mesh_dt = %f \n", 
                  // time, pmb->pmy_mesh->time, dt, pmb->pmy_mesh->dt);
              //   }
              // } 
            }
          }

        }

        // 进入黑洞的处理
        //? 这里的处理方式合适吗？动量为 0 应该没错，能量（压强）应该是个小值还是 0 ？
        if (r < R_in) {
          cons(IDN,k,j,i) = std::min(rho_in_BH, rho);
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          // 修改能量
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) = 0.0 + 0.0; 
          }
        }

        //? 外部区域的处理，应该怎么做？目前是设定为保持初始值
        // TODO 或许可以改为继承该格子之前的值？但怎么操作呢？在 source term 里面不知道能不能做这件事
        // 应该可以把 ProblemGenerator 里的函数抽象出来，然后在这里调用
        if (r > R_out) {
          cons(IDN,k,j,i) = rho_init;
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
  
  // 根据 debug 信息发现，自定义源项在每个时间步会被 call 两次（应该是取决于时间积分器）。
  //* 对于 VL2：第一次是在时间步开始时，传入函数的 dt 为真实 dt 的 1/2（预报步）；另一次在时间步的一半处，传入的是完整的 dt（校正步）。
  // 对于那些 *dt 的增量更新操作不必考虑，但对于直接附加的量（SN 爆炸能量），需要专门处理。
  if (debug >= DEBUG_Mesh && pmb->gid == 0 ){
    printf_and_save_to_stream(debug_stream, 
    "DEBUG_Mesh: Source Term is called at time = %f, mesh_time = %f, dt = %f, mesh_dt = %f \n", 
    time, pmb->pmy_mesh->time, dt, pmb->pmy_mesh->dt);
  }

  return;
}



//TODO 在注入 SuperNova 之后，立即减小 TimeStep？
Real MyTimeStep(MeshBlock *pmb) {
  Real min_dt = std::numeric_limits<Real>::infinity();  // 先初始化为一个很大的数

  if (cooling.cooling_flag) {
    Real cooling_dt = cooling.CoolingTimeStep(pmb);
    min_dt = std::min(min_dt, cooling_dt);
  }
  return min_dt;
}

//========================================================================================
// 自定义 hst Output
//========================================================================================
// hst_NaN 用于在某些配置下，直接输出 nan。这样在读取时就可以很容易发现相应的 hst 输出是不启用的。
// 在 EnrollUserHistoryOutput 时，使用三目运算符来使用 hst_NaN，这样应该比在相应的 hst function 里返回 NaN 要更快一点。
inline Real hst_NaN(MeshBlock *pmb, int iout) {
  return std::numeric_limits<Real>::quiet_NaN();
}

// 统计总共有多少个 MeshBlocks
inline Real hst_num_MeshBlocks(MeshBlock *pmb, int iout) {
  return 1;  // 每个 MeshBlock 返回 1，最后加起来。
}

inline Real hst_dt_hyperbolic(MeshBlock *pmb, int iout) {
  return pmb->pmy_mesh->dt_hyperbolic;
}
inline Real hst_dt_user(MeshBlock *pmb, int iout) {
  return pmb->pmy_mesh->dt_user;  
  // 如果没 EnrollUserTimeStepFunction，这里的 dt_user 会是 Athena++ 内部初始化时设定的 std::numeric_limits<Real>::max()
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


  // 从 input file 中读取参数，存到 pgen 文件的全局变量中
  GM_BH = pin->GetReal("problem","GM_BH");
  R_in = pin->GetReal("problem","R_in");
  R_out = pin->GetReal("problem","R_out");
  rho_in_BH = pin->GetReal("problem","rho_in_BH");
  rho_init = pin->GetReal("problem","rho_init");
  E_tot_init = pin->GetReal("problem", "E_tot_init");

  integrator = pin->GetOrAddString("time", "integrator", "vl2"); // 记录 integrator。目前用来处理超新星的注入时机

  cooling = Cooling(pin, this);

  // 读取超新星的参数
  SN_flag = pin->GetOrAddInteger("supernova","SN_flag",0);
  if (SN_flag > 0) {
    if ( !NON_BAROTROPIC_EOS ) {
      throw std::runtime_error("SN explosion is not implemented for barotropic EOS!"); 
    }
    supernova_list = read_supernova_list(pin, ndim); // 需要传入 ndim，从而确定区域的维数
  }

  // 读取自定义的 AMR 参数
  if (adaptive) {
    SetRootLevel(root_level);
    InitRefinementCondition(pin);
  }

  // <debug> 里的参数
  debug = pin->GetOrAddInteger("debug", "debug", 0);
  if (debug > DEBUG_NONE) {
    debug_filepath = pin->GetOrAddString("debug", "debug_filepath", "info/debug_info.txt");
    ensure_parent_directory_exists(debug_filepath);
    debug_stream.open(debug_filepath, std::ios::app); // 用于输出 debug 信息的文件流
    printf("DEBUG: saving debug info to file %s \n", debug_filepath.c_str());
  }

  // Enroll 各种自定义函数

  // 将自定义的源项注册到 Athena++ 中
  // 调用时机：在每个时间步中调用两次，一次在开头，一次在中间
  EnrollUserExplicitSourceFunction(SMBH_grav);

  // 将自定义的 AMR 条件注册到 Athena++ 中
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  // 自定义的步长限制
  EnrollUserTimeStepFunction(MyTimeStep);  // 无论如何都加入自定义步长限制。如果没开 Cooling 之类的，则 MyTimeStep 立即返回一个很大的数，不影响速度


  // 自定义的 hst Output
  AllocateUserHistoryOutput(3);
  EnrollUserHistoryOutput(0, hst_num_MeshBlocks, "num_MeshBlocks");
  EnrollUserHistoryOutput(1, hst_dt_hyperbolic, "dt_hyperbolic", UserHistoryOperation::min);
  EnrollUserHistoryOutput(2, hst_dt_user, "dt_user", UserHistoryOperation::min); // 有必要把 dt_cooling 单独区分出来的吗？

  // Print unit 相关信息
  punit->PrintCodeUnits();
  punit->PrintConstantsInCodeUnits();

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

  //TEMP heart 的初始条件，等写好了 region 类之后整理
  if (init_cond_type == "heart") {
    std::string block = "initial_condition/heart";
    Real a = pin->GetOrAddReal(block, "a", 3.3);
    Real b = pin->GetOrAddReal(block,"b", 0.75);
    Real heart_size = pin->GetOrAddReal(block,"heart_size", 1.0);
    Real rho_inside = pin->GetReal(block,"rho_inside");
    Real E_tot_inside = pin->GetReal(block,"E_tot_inside");

    Heart heart({0, 0, 0}, heart_size, a, b); // 构造一个心形区域。目前默认 center 是 (0,0,0)

    Real x,y,z;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          z = pcoord->x3v(k);
          y = pcoord->x2v(j);
          x = pcoord->x1v(i);

          if ( heart.contains({x, y, z}) ) {
            phydro->u(IDN,k,j,i) = rho_inside;
            phydro->u(IEN,k,j,i) = E_tot_inside;
          } else {
            phydro->u(IDN,k,j,i) = rho_init;
            phydro->u(IEN,k,j,i) = E_tot_init;
          }

          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;


        }
      }
    }
  };

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


// 调用时机：每个将要输出 output 的时间步的末尾。//? 和 Mesh::UserWorkInLoop 谁先？
// 函数用途：计算用户定义的输出变量
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  if (gid == 0){
    // 自定义的 print 信息： 用 ✅ 标记是否在这一步进行了文件 Output
    std::string spaces = std::string(66, ' '); // 用于跳到 Athena++ 原本的输出的末尾
    std::string marker = "✅" ;

    // 匿名函数，用于在一行内输出 output 的信息
    auto output_in_line = [=](){ // 用 = 表示捕获所有外部变量的值
      std::cout << spaces 
      << marker << " " << pin->GetInteger("output2", "file_number")  // 文件编号。//TODO 目前 hard code 为 output2，以后改为最主要的那个 output
      << " " << pmy_mesh->nbtotal << " " << pmy_mesh->nbnew << " " << pmy_mesh->nbdel   // MeshBlock 的数量
      << "\r"; // 输出完后回到行首，让 Athena++ 打印其原本要打印的信息
    };

    if (pmy_mesh->time == 0.0) { // 第一步因为 Athena++ 会输出 "\n Setup complete ... \n" 什么的，所以要特殊处理
      std::cout << "\n" << spaces << "Output" << std::string(2, ' ') << " MeshBlocks"; // header 跟 Setup complete ... 同一行
      std::cout << std::string(2, '\n');     
      output_in_line();
      std::cout << "\033[A\033[A\033[A";        // 输出完后再向上 3 行回去
    } else {
      output_in_line();
    }
    //TODO 在输出 output 时，也把当前时间输出到 ./info 下一个单独的文件，这样就不用从每个 ds 中单独读取了。但需要找到 pin 中合适的 <output[n]> block
  }
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
