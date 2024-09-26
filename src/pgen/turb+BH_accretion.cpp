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
#include "../scalars/scalars.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

/* -------------------------------------------------------------------------- */
/*                            自定义的变量、函数等                               */
/* -------------------------------------------------------------------------- */

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
#include "ckj_code/my_outputs.hpp"


// 用于在 SourceTerm 中判断当前是否为最后一个 stage。目前使用 dt_ratio 来判断。
inline bool IsLastStage(const Real dt_ratio) {
  // 目前看来，我猜测 orbital_advection 和 integrator 不能在运行时动态修改，但可以在 restart 时通过 cmd 或 input file 更改。因此，只需在 InitUserMeshData 中拿到 integrator 即可。但有待验证。
  return dt_ratio == beta_last_stage;
}

// 这个函数只在 time integrator 的最后一个 stage 调用。通过 IsLastStage 函数来判断。
void SourceTermAtLastStage(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  //* 在这个函数中，传给下层 SourceTerm 的 time 和 dt 应该是完整的 mesh_time 和 mesh_dt，而非 SourceTerm 中传入的当前 stage 的 dt。
  // 例如 Cooling，需要在这里走一整步。而 SN SourceTerm 内部的 dt/mesh_dt 也得到 1.0，意味着在这里一次性注入。
  // 在下层的 SourceTerm 中应该不需要用到「当前 stage」的 time 和 dt，相关的判断都在外层做完了。
  const Real mesh_time = pmb->pmy_mesh->time;
  const Real mesh_dt = pmb->pmy_mesh->dt;

  // 把 SN 放在 cooling 之后，避免 cooling 时受到 SN 的影响，因为这一步尚未根据 SN 限制 dt。

  // Cooling，在这里意味着使用当前步开始时的 prim，只走一整步。//? 这可能不一定稳定 or 准确？
  if (cooling.cooling_flag && cooling.source_term_position == SourceTermPosition::AfterSourceTerm) {
    cooling.CoolingSourceTerm(pmb, mesh_time, mesh_dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // 超新星的注入 //*目前默认是在 SourceTerm 的结尾。
  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::AfterSourceTerm) {
    supernovae.SuperNovaeSourceTerm(pmb, mesh_time, mesh_dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }
};

// Source Term 源项
//* 应该只依赖于 prim 而修改 cons。不能修改 prim。不应该依赖于 cons（尤其是因为 cons 可能被前面的源项修改）。
// 把各个源项拆分。这样需要反复进入循环好几次，但是更加清晰。也可以考虑把最内层循环内的各个步骤抽象成函数？
void MySourceFunction(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  // CoolingSourceTerm 与其它 SourceTerm 的顺序先后应该并无影响，因为反正都是根据这一步未修改的 prim 去更新 cons。而且 CoolingTimeStep 也是在这一步开始前就计算好的。
  // 但是，第一个 stage 对 cons 的更新会进入第二个 stage 的 prim，
  if (cooling.cooling_flag && cooling.source_term_position == SourceTermPosition::InSourceTerm) {
    cooling.CoolingSourceTerm(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // 超新星的注入 // 目前默认并不在这里注入，如果随着源项注入，可能会炸步长 / 注入比例不对，但未经测试（暂时没必要尝试这个）。
  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::InSourceTerm) {
    supernovae.SuperNovaeSourceTerm(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // 中心黑洞的引力
  if (M_BH != 0.0) {
    SMBH_grav(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  Real mesh_dt = pmb->pmy_mesh->dt;
  if (IsLastStage(dt/mesh_dt)) { // 若为最后一个 stage，则在这里插入一些操作：比如 SN 的注入、Cooling 的操作等
    SourceTermAtLastStage(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

}

void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
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
        vx  = prim(IVX,k,j,i);
        vy  = prim(IVY,k,j,i);
        vz  = prim(IVZ,k,j,i);

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
        }

        //TODO 把是否在黑洞内、是否在 R_out 之外的判断抽出来，放到一个 SimulationRegion 类的实例方法中？还是放到 SinkRegion 类中？
        //TODO 把黑洞内的处理、外部区域的处理（保持初值）抽出来，放到 SimulationRegion 类中。这两个函数的调用是不是应该放在 AfterSourceTerm 中？这样只需调用一次，而且 RK2 关心的是变化量，可能不能准确一步抵达目标？
        //* 如果要把 sink 内的处理抽出来，要注意源项的顺序，因为 sink 的源项要记录吸积物理量，依赖于 cons
        // 进入黑洞的处理
        if (r < R_in) {
          // 记录被 sink region 吸收的物理量：质量、角动量、SN tracer。
          //* 注意这里不要用 prim，而是用 cons，因为 prim 尚未更新，是上一步结尾的值！而且要先记录这些量的改变量，再修改 cons。
          const Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);    //FUTURE 改用 pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol); ？ 这样可以利用 OpenMP SIMD 一次算一行。但如果不是性能热点，可能没必要。
          const Real cell_mass = cons(IDN,k,j,i) * cell_volume;

          // 质量
          //BUG 这里把 density floor 的质量也一并算进去了，最好是剔除。不过影响不太大
          pmb->ruser_meshblock_data[I_accreted_mass](0) += cell_mass;  // 记录被 BH 吸收的质量
          pmb->ruser_meshblock_data[I_accretion_rate](0) += cell_mass;  // 同时也记到 I_accretion_rate 里，最后加总后再除以 dt。注意这里必须是 +=，因为这是对每个格点的质量进行累加。

          // 动量
          Vector momentum;
          for (int l = 0; l < 3; ++l) {
            momentum[l] = cons(IM1+l,k,j,i);
            pmb->ruser_meshblock_data[I_accreted_momentum](l) += momentum[l];  // 记录被 BH 吸收的动量
          }

          // 角动量
          Vector angular_momentum = CrossProduct({x,y,z}, momentum);   //* 目前假设 BH 位于 {0,0,0}，否则 {x,y,z} 要减去 BH 的位置，而且速度也要改为相对速度！
          for (int l = 0; l < 3; ++l) {
            pmb->ruser_meshblock_data[I_accreted_angular_momentum](l) += angular_momentum[l];  // 记录被 BH 吸收的角动量
          }
          
          // SN Tracer
          if (NSCALARS > 0) {
            Real SN_tracer_mass = cons_scalar(PassiveScalarIndex::SN,k,j,i) * cell_volume;
            pmb->ruser_meshblock_data[I_accreted_SN_tracer](0) += SN_tracer_mass;
            pmb->ruser_meshblock_data[I_accretion_rate_SN_tracer](0) += SN_tracer_mass;
            cons_scalar(PassiveScalarIndex::SN,k,j,i) = 0;  // 清空 SN Tracer
          }

          // 处理 Sink Region 内的物理量改变
          //? 这里似乎不需要 min，因为 sink 内的流体密度（= sink密度 + 新进来的流体的贡献）总是高于我所设定的 sink 密度。但暂时先保留。如果观察到 sink 内密度降得比设定值还低，也许有问题。
          cons(IDN,k,j,i) = std::min(rho_in_BH, rho);  // 使用 min，如果密度低于设定值，则仍保持其密度。但按说不会发生。
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          // 修改能量
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) = 0.0 + 0.0;  // 把内能（压强）清零（压强会由于 floor 而有一个小值）。动能也是 0.
          }
        }
        //? 外部区域的处理，应该怎么做？目前是设定为保持初始值
        //? 或者，应该用边界处理？看看其他论文怎么做的，问问 Kohei。
        //* 也许还需要允许选择不同的处理方式？不一定非得保持初值？尤其是对于 approx_Bondi 这种。
        // TODO 或许可以改为继承该格子之前的值？但怎么操作呢？在 source term 里面不知道能不能做这件事
        // 应该可以把 ProblemGenerator 里的函数抽象出来，然后在这里调用
        if (r > R_out) {
          initial_condition->SetSingleCell(pmb, i, j, k);
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
// 自定义的 dt_user。注意这个函数的返回值应该严格 > 0 ，否则虽然 Integrator 不报错，但感觉很危险…
Real MyTimeStep(MeshBlock *pmb) {
  Real min_dt = std::numeric_limits<Real>::max();  // 先初始化为一个很大的数。这里不用 infinity() 是为了万一启用 -ffast-math 之类的编译选项的话，对无穷大的处理可能不正确
  //* 如果 SN 在 UserWorkInLoop 中注入，那么由于在决定下一个 cycle 的 dt 时尚未注入 SN，因此无法自动限制步长，需要手动限制：
  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::UserWorkInLoop) {
    // 根据 supernova_to_inject 是否为空来判断，实际上会把所有 MeshBlock 的都算进来。但反正是取最小值，所以不影响结果。
    if (!supernovae.supernova_to_inject.empty()) {
      min_dt = std::min(min_dt, 1e-8);  //TEMP 这里的数值是手动写的，可能需要更改
    }
  }

  if (cooling.cooling_flag) {
    Real cooling_dt = cooling.CoolingTimeStep(pmb);  // 注：这里使用的是 cons2prim 之后的 prim 来计算 cooling_dt
    min_dt = std::min(min_dt, cooling_dt);
  }
  return min_dt;
}


/* -------------------------------------------------------------------------- */
/*                            自定义 hst Output                                */
/* -------------------------------------------------------------------------- */


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

inline Real hst_accreted_mass(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_mass](0);
}
// 输出当前时间步内的 accretion rate
//BUG hst enroll 的函数每 dcycle 步才调用一次！不应该在这里清零！这造成了 accretion rate 的输出大了 dcycle 倍！
//* 修这个 bug 时，要注意 hst 函数间隔 dcycle 才调用、各节点独立、 InSourceTerm 被调用两次且 dt 不同、累加清空时机
inline Real hst_accretion_rate(MeshBlock *pmb, int iout) {
  Real accretion_rate = pmb->ruser_meshblock_data[I_accretion_rate](0) / pmb->pmy_mesh->dt;  // 除以 dt，得到每秒的 accretion rate
  pmb->ruser_meshblock_data[I_accretion_rate](0) = 0;  // 重置
  return accretion_rate;
}

inline Real hst_accreted_SN_tracer(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_SN_tracer](0);
}
inline Real hst_accretion_rate_SN_tracer(MeshBlock *pmb, int iout) {
  Real accretion_rate_SN_tracer = pmb->ruser_meshblock_data[I_accretion_rate_SN_tracer](0) / pmb->pmy_mesh->dt;  // 除以 dt，得到每秒的 accretion rate
  pmb->ruser_meshblock_data[I_accretion_rate_SN_tracer](0) = 0;  // 重置
  return accretion_rate_SN_tracer;
}


// 动量的三个分量
inline Real hst_accreted_momentum_x(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_momentum](0);
}
inline Real hst_accreted_momentum_y(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_momentum](1);
}
inline Real hst_accreted_momentum_z(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_momentum](2);
}

// 角动量的三个分量
inline Real hst_accreted_angular_momentum_x(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_angular_momentum](0);
}
inline Real hst_accreted_angular_momentum_y(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_angular_momentum](1);
}
inline Real hst_accreted_angular_momentum_z(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[I_accreted_angular_momentum](2);
}



/* -------------------------------------------------------------------------- */
/*                       以下是 Athena++ 提供的接口                              */
/* -------------------------------------------------------------------------- */


// 调用时机：Mesh 类实例化时，也即 main.cpp 的 Step 4。restart 时也会调用。
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
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << '\n'
        << "non zero Turbulence flag is set without FFT!" << '\n';
    ATHENA_ERROR(msg);
    return;
#endif
  }

  // 用于 IsLastStage 函数
  if (orbital_advection != 0) {
    std::cout << "### ERROR: 启用了 orbital_advection，但目前我对 SourceTermAtLastStage 的实现只使用 dt_ratio 的值来 \
    简单判断。对于 orbital_advection，需要拿到 TaskList 的指针，获取最后一个 stage 的 beta；或者直接根据 stage 判断。\n";
    throw std::runtime_error("orbital_advection 与自定义的 SourceTermAtLastStage 不兼容！");
  }
  // 以下各个 Integrator 的 beta 的值都是假设没有启用 orbital_advection 的。从 time_integrator.cpp 中抄来：
  // rk1: beta == {1.0}
  // vl2: beta == {0.5, 1.0}
  // rk2: beta == {1.0, 0.5}
  // rk3: beta == {1.0, 0.25, TWO_3RD}
  // rk4: beta == {1.193743905974738, 0.099279895495783, 1.131678018054042, 0.310665766509336}
  // ssprk5_4: 太长了……反正现在也用不着
  std::string integrator = pin->GetString("time", "integrator");
  std::unordered_map<std::string, Real> beta_last_stage_map = {
    {"rk1", 1.0},
    {"vl2", 1.0},
    {"rk2", 0.5},
    {"rk3", TWO_3RD},
  };
  beta_last_stage = beta_last_stage_map.at(integrator);  // 使用 .at 来获取值，如果 key 不存在会抛出异常。

  // 从 input file 中读取参数，存到 pgen 文件的模块级别变量中
  Real M_BH_in_Msun = pin->GetOrAddReal("problem","M_BH", 0.0);
  M_BH = M_BH_in_Msun * punit->solar_mass_code;
  GM_BH = M_BH * punit->grav_const_code;

  R_in = pin->GetOrAddReal("problem","R_in", 0.0);
  R_out = pin->GetOrAddReal("problem","R_out", std::numeric_limits<Real>::max());
  rho_in_BH = pin->GetReal("problem","rho_in_BH");  //TODO 把 rho_in_BH 改名为 n_in_BH，设定默认值，并且在默认模板中不设定值。

  // 构造我自定义机制的对象
  //* 注意：在这里构造的对象，都必须是确定性的，从而保证在 restart 时 / 跨 MPI rank 的一致性。
  initial_condition = std::unique_ptr<InitialCondition>(new InitialCondition(this, pin));

  cooling = Cooling(this, pin);  // 初始化 cooling

  supernovae = Supernovae(this, pin);  // 初始化 supernovae

  if (supernovae.source_term_position <= SourceTermPosition::AfterSourceTerm && cooling.source_term_position == SourceTermPosition::UserWorkInLoop) {
    std::cout << "### WARNING: Supernova 在 SourceTerm 中/结尾 注入，而 Cooling 在 UserWorkInLoop 中。这会导致注入 SN 那一步的 Cooling TimeStep 未受限制！ \n";
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
  // 调用时机：根据 integrator 有所不同。对于 VL2，在每个时间步中调用两次，一次在开头，一次在中间；
  EnrollUserExplicitSourceFunction(MySourceFunction);

  // 将自定义的 AMR 条件注册到 Athena++ 中
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  // 自定义的步长限制
  EnrollUserTimeStepFunction(MyTimeStep);  // 无论如何都加入自定义步长限制。如果没开 Cooling 之类的，则 MyTimeStep 立即返回一个很大的数，不影响速度


  // 自定义的 hst Output
  enum UserHistoryOutputIndex {
    I_num_MeshBlocks,
    I_dt_hyperbolic,
    I_dt_user,
    I_accreted_mass,
    I_accretion_rate,
    I_accreted_SN_tracer,
    I_accretion_rate_SN_tracer,
    I_accreted_momentum_x,
    I_accreted_momentum_y,
    I_accreted_momentum_z,
    I_accreted_angular_momentum_x,
    I_accreted_angular_momentum_y,
    I_accreted_angular_momentum_z,
    // 其他 ...
    N_UserHistoryOutput // 总个数
  };

  // 如果不 Enroll 某些 hst 函数，那么对应的 user_history_func_[i] 就会是 nullptr，而在 HistoryOutput 中有判空。 //? 会输出 nan 还是 0 还是什么？有待试验
  //TODO 给这些 hst 加条件，比如 accretion 相关量只有在存在 sink 的时候启用；SN 相关量只有在存在 SN 的时候启用。
  //? 用 hst_NaN 来给出空输出，还是直接用默认的 nullptr 的输出？

  AllocateUserHistoryOutput(N_UserHistoryOutput);
  EnrollUserHistoryOutput(I_num_MeshBlocks, hst_num_MeshBlocks, "num_MeshBlocks");
  EnrollUserHistoryOutput(I_dt_hyperbolic, hst_dt_hyperbolic, "dt_hyperbolic", UserHistoryOperation::min);
  EnrollUserHistoryOutput(I_dt_user, hst_dt_user, "dt_user", UserHistoryOperation::min); // 有必要把 dt_cooling 单独区分出来的吗？
  EnrollUserHistoryOutput(I_accreted_mass, hst_accreted_mass, "accreted_mass");
  EnrollUserHistoryOutput(I_accretion_rate, hst_accretion_rate, "accretion_rate");
  EnrollUserHistoryOutput(I_accreted_SN_tracer, hst_accreted_SN_tracer, "accreted_SN_tracer");
  EnrollUserHistoryOutput(I_accretion_rate_SN_tracer, hst_accretion_rate_SN_tracer, "accretion_rate_SN_tracer");
  EnrollUserHistoryOutput(I_accreted_momentum_x, hst_accreted_momentum_x, "accreted_momentum_x");
  EnrollUserHistoryOutput(I_accreted_momentum_y, hst_accreted_momentum_y, "accreted_momentum_y");
  EnrollUserHistoryOutput(I_accreted_momentum_z, hst_accreted_momentum_z, "accreted_momentum_z");
  EnrollUserHistoryOutput(I_accreted_angular_momentum_x, hst_accreted_angular_momentum_x, "accreted_angular_momentum_x");  
  EnrollUserHistoryOutput(I_accreted_angular_momentum_y, hst_accreted_angular_momentum_y, "accreted_angular_momentum_y");
  EnrollUserHistoryOutput(I_accreted_angular_momentum_z, hst_accreted_angular_momentum_z, "accreted_angular_momentum_z");

  if (Globals::my_rank == 0) {
    std::cout << "pgen 编译于: " << __DATE__ << " " << __TIME__ << '\n'; 
    // Print unit 相关信息
    punit->PrintCodeUnits();
    punit->PrintConstantsInCodeUnits();
    //TODO Print 自定义的参数
  }

  return;
}

// 调用时机：程序启动时（包括 restart 时）
// 函数用途：初始化每个 MeshBlock 上的局部变量（包括数组数据）
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // 之前我用来传递 input 参数，后来改成了用 pgen 的全局变量
  AllocateRealUserMeshBlockDataField(N_RealUserMeshBlockData);
  ruser_meshblock_data[I_accreted_mass].NewAthenaArray(1);  // 黑洞吸积的总质量，在每个 MeshBlock 上单独存储。
  ruser_meshblock_data[I_accretion_rate].NewAthenaArray(1);  // 用于计算黑洞吸积率，储存当前 dt 内吸积的质量，最后除以 dt。
  ruser_meshblock_data[I_accreted_SN_tracer].NewAthenaArray(1);  // 被黑洞吸积的 SN ejecta 的质量，在每个 MeshBlock 上单独存储。
  ruser_meshblock_data[I_accretion_rate_SN_tracer].NewAthenaArray(1);  // 用于计算吸积率中 SN tracer 的质量。储存当前 dt 内吸积的质量，最后除以 dt。
  ruser_meshblock_data[I_accreted_momentum].NewAthenaArray(3);  // 黑洞吸积的总动量，在每个 MeshBlock 上单独存储。
  ruser_meshblock_data[I_accreted_angular_momentum].NewAthenaArray(3);  // 黑洞吸积的总角动量，在每个 MeshBlock 上单独存储。

  // 设定 uov
  AllocateUserOutputVariables(N_UOV);  // allocate 我自定义的 output 变量（uov, user_out_var）
  //* 为每个 user_out_var 变量设置名字。这里使用 _snake_case，额外在开头加一个 _ 以规避 yt 中名称重复导致的不便。
  SetUserOutputVariableName(I_cooling_rate, "_cooling_rate");


  return;
}


// 调用时机：程序启动时（应该不包括 restart 时？）
// 函数用途：设定初始条件（每个 MeshBlock 的）
void MeshBlock::ProblemGenerator(ParameterInput *pin) {  
  initial_condition->SetInitialCondition(this);
}


// 调用时机：在每个时间步的结尾。每个 MeshBlock 调用一次
// 函数用途：仅用于分析，不用于操作数据
void MeshBlock::UserWorkInLoop() {
  //* 如果这里以后要实现多个源项的作用，要考虑其顺序：
  //? 是像 SourceTerm 里那样全都从同一个 prim 出发，改变 cons，最后再统一 cons2prim； 还是进行 operator splitting，一项处理完之后先 cons2prim 一下，再进行下一项？目前采用前者，都从同一个 prim 出发。

  if (cooling.cooling_flag && cooling.source_term_position == SourceTermPosition::UserWorkInLoop) {  // 如果开启了 operator splitting，则在 UserWorkInLoop 中调用 CoolingSourceTerm
    // 这里调用的参数，参考 src/task_list/time_integrator.cpp 中 AddSourceTerms 的调用
    cooling.CoolingSourceTerm(this, pmy_mesh->time, pmy_mesh->dt,
                              phydro->w, pscalars->r, pfield->bcc,
                              phydro->u, pscalars->s);
  }

  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::UserWorkInLoop) {
    supernovae.SuperNovaeSourceTerm(this, pmy_mesh->time, pmy_mesh->dt,
                                  phydro->w, pscalars->r, pfield->bcc,
                                  phydro->u, pscalars->s);
  }

  // 在这里做 cons2prim 的话，由于 MeshBlock 之间的通信在此之前，会造成数据不自洽，导致 MeshBlock 的边界上不对。因此这里不能 cons2prim

  // 目前还不清楚，但我看其他的几个 pgen 里也有在这里自己施加 floor 的
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


// 调用时机：每个将要输出 output 的时间步的末尾。(但不包括 hst 输出，因为太频繁了) //? 和 Mesh::UserWorkInLoop 谁先？
// 函数用途：计算用户定义的输出变量
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // 计算 uov
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        // user_out_var(I_test,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);  // 例子
        // cooling_rate 不在这里计算，而是每一个 timestep 在计算 cooling 时储存到 uov 中。
      }
    }
  }


  // 自定义一些实时输出
  if (gid == 0){
    // 自定义的 print 信息： 用 ✅ 标记是否在这一步进行了文件 Output
    std::string spaces = std::string(66, ' '); // 用于跳到 Athena++ 原本的输出的末尾
    std::string marker = "✅" ;

    // 匿名函数，用于在一行内输出 output 的信息
    auto output_in_line = [&](){ // 用 = 表示捕获所有外部变量的值
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
  }
  return;
}


// 调用时机：模拟结束时
// 函数用途：清理 MeshBlock::InitUserMeshData 分配的资源
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // 这个函数用于在模拟结束后做一些事情
  if (debug >= DEBUG_Main) {
    printf_and_save_to_stream(debug_stream, "DEBUG_Main: Mesh::UserWorkAfterLoop is called \n");
  }

  // 关闭 debug_stream
  if (debug_stream.is_open()) { 
    printf("debug_stream is open, closing \n");
    debug_stream.close();
  }
}
