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
#include <chrono>
#include <iomanip>

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

// 自定义的头文件（各种模块）
#include "ckj_code/ckj_plugin.hpp"
#include "ckj_code/utils.hpp"
#include "ckj_code/region.hpp"
#include "ckj_code/refinement_condition.hpp"
#include "ckj_code/supernova.hpp"
#include "ckj_code/cooling.hpp"
#include "ckj_code/my_outputs.hpp"
#include "ckj_code/hst_output.hpp"
#include "ckj_code/progress_report.hpp"

// 本 pgen 的头文件（以后考虑合并）
#include "turb+BH_accretion.hpp"


/* -------------------------------------------------------------------------- */
/*                                Source Terms                                */
/* -------------------------------------------------------------------------- */

// 把各个源项拆分。这样需要反复进入循环好几次，但是更加清晰。
void MySourceFunction(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  // 在自定义源项的开头 CheckCell
  if (debug >= DEBUG_Cell) {
    CheckCell(false, pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

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

  if (ckj_plugin::current_stage == ckj_plugin::nstages) { // 若为最后一个 stage，则在这里插入一些操作：比如 SN 的注入、Cooling 的操作等
    SourceTermAtLastStage(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // SMBH Sink Region 的处理，需要放在最后，因为记录吸积量依赖于 cons。
  if (R_in > 0.0) {
    SMBH_sink(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // 在自定义源项的结尾 CheckCell
  if (debug >= DEBUG_Cell) {
    CheckCell(true, pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

}


void SourceTermAtLastStage(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  //* 在这个函数中，传给下层 SourceTerm 的 time 和 dt 应该是完整的 mesh_time 和 mesh_dt，而非 SourceTerm 中传入的当前 stage 的 dt。
  // 例如 Cooling，需要在这里走一整步。而 SN SourceTerm 内部的 dt/mesh_dt 也得到 1.0，意味着在这里一次性注入。
  // 在下层的 SourceTerm 中应该不需要用到「当前 stage」的 time 和 dt，相关的判断都在外层做完了。
  const Real mesh_time = pmb->pmy_mesh->time;
  const Real mesh_dt = pmb->pmy_mesh->dt;

  // 把 SN 放在 cooling 之后，避免 cooling 时受到 SN 的影响，因为这一步尚未根据 SN 限制 dt。虽然理论上源项不应该依赖于 cons，但 cooling 可能依赖于 cons 来判断还剩多少能量？

  // Cooling，在这里意味着使用当前步开始时的 prim，只走一整步。//? 这可能不一定稳定 or 准确？
  if (cooling.cooling_flag && cooling.source_term_position == SourceTermPosition::AfterSourceTerm) {
    cooling.CoolingSourceTerm(pmb, mesh_time, mesh_dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  // 超新星的注入 //*目前默认是在 SourceTerm 的结尾。
  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::AfterSourceTerm) {
    supernovae.SuperNovaeSourceTerm(pmb, mesh_time, mesh_dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

};


void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  
  const Real source_weight = ckj_plugin::source_term_weight; // 在循环外把这个值存下来，避免每次循环都进行命名空间、全局变量的查找
  Real BH_gravity_work = 0.0;  // 用于在 SIMD 循环中累计

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);
#pragma omp simd reduction(+:BH_gravity_work)  // SIMD 矢量化。必须要加 reduction 语句，否则 BH_gravity_work 会被多个线程同时写入，导致结果偏小！
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real r = std::sqrt(x*x + y*y + z*z);
        Real r3 = r*r*r;
        Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);
        
        Real rho = prim(IDN,k,j,i); 
        Real vx  = prim(IVX,k,j,i);
        Real vy  = prim(IVY,k,j,i);
        Real vz  = prim(IVZ,k,j,i);

        // 只对黑洞外的区域施加引力。在这里对 R_in 以内施加引力没用，接下来会被重设
        //? 这里需要 && r <= R_out 吗？
        if (R_in <= r and r <= R_out) {  // 这里限定为 R_in 和 R_out 之间，因为其他地方会被覆盖，引力做功的计算就不准确了
          cons(IM1,k,j,i) += - GM_BH*x/r3 * rho * dt;
          cons(IM2,k,j,i) += - GM_BH*y/r3 * rho * dt;
          cons(IM3,k,j,i) += - GM_BH*z/r3 * rho * dt;

          // 能量的改变：非全局正压时需要更新能量密度
          if (NON_BAROTROPIC_EOS) {
            // 这里可以稍微优化，快一点
            Real dE = - GM_BH*(x*vx + y*vy + z*vz)/r3 * rho * dt;
            cons(IEN,k,j,i) += dE; // 根据动量的改变，相应地改变动能。内能目前不变。
            //? 这里是否需要考虑把平方项补偿上去？有待测试。
            BH_gravity_work += dE * cell_volume * source_weight;  // 记录引力做的功
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

  auto hst_data = hst->get_proxy(pmb);
  hst_data["BH_gravity_work"] = BH_gravity_work;  // 用于输出到 hst
  
  return;
}



//FUTURE 放到 SimulationRegion 类中？
void SMBH_sink(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  auto hst_data = hst->get_proxy(pmb);

  // 在第一个 stage 清零吸积率，用于统计当前 cycle 内的吸积率。目前放在 SMBH_sink 源项内，是为了把涉及吸积量的计算都放在一起。
  if (ckj_plugin::current_stage == 1) { 
    hst_data["accretion_rate"] = 0.0;
    hst_data["accretion_rate_SN_tracer"] = 0.0;
  }

  Real mesh_dt = pmb->pmy_mesh->dt;
  const Real source_weight = ckj_plugin::source_term_weight; // 在循环外把这个值存下来，避免每次循环都进行命名空间、全局变量的查找
  
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) { // 由于这个函数涉及许多吸积量的记录，是 reduction 操作，所以不启用 SIMD 矢量化。
        Real x = pmb->pcoord->x1v(i);
        Real r = std::sqrt(x*x + y*y + z*z);

        Real dens = cons(IDN,k,j,i); // 用 dens 称呼 cons 守恒密度，与 prim 密度区分。

        //FUTURE 把是否在黑洞内、是否在 R_out 之外的判断抽出来，放到一个 SimulationRegion 类的实例方法中？还是放到 SinkRegion 类中？
        // 进入黑洞的处理
        if (r < R_in) {
          // 记录被 sink region 吸收的物理量：质量、角动量、SN tracer。
          //* 注意这里不要用 prim，而是用 cons，因为 prim 尚未更新，是上一步结尾的值！而且要先记录这些量的改变量，再修改 cons。
          const Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);    //FUTURE 改用 pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol); ？ 这样可以利用 OpenMP SIMD 一次算一行。但如果不是性能热点，可能没必要。其他 pgen 里用的 vol 长度都是 pmb->ncells1，但这个是包含了 2*NGHOST 的，感觉 pmb->blocksize->nx1 就够了，不知道为什么。
          const Real dens_new = std::min(rho_sink, dens);  // 使用 min，如果密度低于 rho_sink，则仍保持其密度，而非抬升到 rho_sink。
          //* 考虑到有强风（比如 SN shock）吹过时，sink 周围的径向速度可能是正的。尽管 sink 内速度为 0，但由于边界值重建，可能产生向外的 flux，导致 sink 内密度减小！
          // 不过，由于 rho_sink 很低，所以即便有强风，也不会向外漏出多少质量。而且一旦吸积，dens 比 rho_sink 低的部分就会很快被填平。

          // 质量
          const Real accreted_mass = (dens - dens_new) * cell_volume * source_weight;       // 计算当前 cell 内减去的质量。根据前面的 min，目前这里是非负的
          hst_data["accreted_mass"] += accreted_mass;             // 记录被 BH 吸收的质量
          hst_data["accretion_rate"] += accreted_mass / mesh_dt;  // 记录吸积率。注意这里分母是 mesh_dt，而不是 dt，因为各个 stage 中传入的 dt 不同，但都应该算是在这个 cycle 对应的 dt 中吸积的质量。

          // 能量
          const Real accreted_energy = cons(IEN,k,j,i) * cell_volume * source_weight;  // 计算当前 cell 内的能量（包括动能和内能）
          hst_data["accreted_energy"] += accreted_energy;  // 记录被 BH 吸收的能量

          // 动量
          Vector momentum;
          for (int l = 0; l < 3; ++l) {
            momentum[l] = cons(IM1+l,k,j,i) * cell_volume * source_weight;
            hst_data[std::string("accreted_momentum_") + "xyz"[l]] += momentum[l]; // 记录被 BH 吸收的动量。注意这里必须用 string 否则两个都是 char，+ 不是字符串连接！
          }

          // 角动量
          Vector angular_momentum = CrossProduct({x,y,z}, momentum);   //* 目前假设 BH 位于 {0,0,0}，否则 {x,y,z} 要减去 BH 的位置，而且速度也要改为相对速度！
          for (int l = 0; l < 3; ++l) {
            hst_data[std::string("accreted_angular_momentum_") + "xyz"[l]] += angular_momentum[l]; // 记录被 BH 吸收的角动量
          }
          
          // SN Tracer
          if (NSCALARS > 0) {
            Real SN_tracer_mass = cons_scalar(PassiveScalarIndex::SN,k,j,i) * cell_volume * source_weight;
            hst_data["accreted_SN_tracer"] += SN_tracer_mass;
            hst_data["accretion_rate_SN_tracer"] += SN_tracer_mass / mesh_dt;  // 记录 SN Tracer 的吸积率
          }


          // 处理 Sink Region 内的物理量改变
          cons(IDN,k,j,i) = dens_new;  // 把密度设定为 sink 密度
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          // 修改能量
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) = 0.0 + 0.0;  // 把内能（压强）清零（压强会由于 floor 而有一个小值）。动能也是 0.
          }
          if (NSCALARS > 0) {
            // 清空所有的 passive scalar，否则会从 sink 中漏出来。推测是由于边界值重建导致 flux 不完全为 0
            // 目前没有针对某些单独的 passive scalar 处理，以后可能添加。
            for (int n=0; n<NSCALARS; ++n) {
              cons_scalar(n,k,j,i) = 0.0;
            }
          }
        }
      }
    }
  }
}


void CheckCell(bool after_source, MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real r = std::sqrt(x*x + y*y + z*z);

        for (std::string check: {"nan", "negative"}) {  // 遍历检查项
          for (std::string vars_str: {"cons", "prim"}) {
            for (int n=0; n <= 4; ++n){
              auto& vars = vars_str == "cons" ? cons : prim;
              std::string position = after_source ? "after" : "before";
              Real value = vars(n,k,j,i);
              
              if (check == "nan" and std::isnan(value) 
              or (check == "negative" and (n == IDN or n == IPR) and value < 0.0)) {
                debug_stream << position << ": " << check << " found in " << vars_str << "("<< n << "," << k << "," << j << "," << i << "), "
                              << "value = " << value << ", ncycle = " << pmb->pmy_mesh->ncycle << ", time = " << time << ", dt = " << dt << ", stage = " << ckj_plugin::current_stage
                               << ", (x,y,z,r) = (" << x << "," << y << "," << z << "," << r << ")\n";
              }
            }
          }
        }
        
      }
    }
  }
}


// 自定义的 dt_user。注意这个函数的返回值应该严格 > 0 ，否则虽然 Integrator 不报错，但感觉很危险…
Real MyTimeStep(MeshBlock *pmb) {
  Real min_dt = std::numeric_limits<Real>::max();  // 先初始化为一个很大的数。这里不用 infinity() 是为了万一启用 -ffast-math 之类的编译选项的话，对无穷大的处理可能不正确
  //* 如果 SN 在 UserWorkInLoop 中注入，那么由于在决定下一个 cycle 的 dt 时尚未注入 SN，因此无法自动限制步长，需要手动限制：
  if (supernovae.SN_flag > 0 && supernovae.source_term_position == SourceTermPosition::UserWorkInLoop) {
    // 根据 supernova_to_inject 是否为空来判断，实际上会把所有 MeshBlock 的都算进来。但反正是取最小值，所以不影响结果。
    if (!supernovae.supernova_to_inject.empty()) {
      min_dt = std::min(min_dt, 1e-8);  // 这里的数值是手动写的，可能需要更改。但反正目前决定 SN 不在 UserWorkInLoop 中注入，因此这里不会启用
    }
  }

  if (cooling.cooling_flag) {
    Real cooling_dt = cooling.CoolingTimeStep(pmb);  // 注：这里使用的是 cons2prim 之后的 prim 来计算 cooling_dt
    min_dt = std::min(min_dt, cooling_dt);
  }
  return min_dt;
}




// 注册所有的 HST 输出条目
void HSTManager::register_all_entries() {
  add_entry("num_MeshBlocks");
  add_entry("dt_hyperbolic", UserHistoryOperation::min);
  add_entry("dt_user", UserHistoryOperation::min);

  //TODO 改为用循环批量添加

  // 黑洞吸积量
  if (R_in > 0.0) {    // 只有当存在 sink region 时启用
    add_entry("accreted_mass");                // 黑洞吸积的总质量
    add_entry("accretion_rate");               // 黑洞当前 cycle 的吸积率。
    add_entry("accreted_SN_tracer");           // 黑洞吸积的 SN ejecta 的质量
    add_entry("accretion_rate_SN_tracer");     // 黑洞当前 cycle 的 SN Tracer 的吸积率。
    add_entry("accreted_energy");              // 黑洞吸积的总能量
    add_entry("accreted_momentum_x");          // 黑洞吸积的总动量
    add_entry("accreted_momentum_y");       
    add_entry("accreted_momentum_z");
    add_entry("accreted_angular_momentum_x");  // 黑洞吸积的总角动量
    add_entry("accreted_angular_momentum_y");
    add_entry("accreted_angular_momentum_z");
  }

  if (cooling.cooling_flag) {
    add_entry("total_cooling_loss");
  }
  if (M_BH != 0.0 and NON_BAROTROPIC_EOS) {
    add_entry("BH_gravity_work");   // SMBH 引力对流体做的功
  }
  if (supernovae.SN_flag > 0) {
    add_entry("SN_injected_energy");  // SN 注入的总能量
    add_entry("SN_injected_mass");    // SN 注入的总质量
    add_entry("SN_injected_number", UserHistoryOperation::min);  // SN 注入的总个数
  }
  // 总角动量
  add_entry("total_angular_momentum_x");
  add_entry("total_angular_momentum_y");
  add_entry("total_angular_momentum_z");

  // multi-phase ISM 的 hst output
  add_entry("ISM_hot_volume");
  add_entry("ISM_warm_volume");
  add_entry("ISM_cold_volume");

  add_entry("ISM_hot_mass");
  add_entry("ISM_warm_mass");
  add_entry("ISM_cold_mass");

  add_entry("ISM_hot_thermal_energy");
  add_entry("ISM_warm_thermal_energy");
  add_entry("ISM_cold_thermal_energy");

  add_entry("ISM_hot_kinetic_energy");
  add_entry("ISM_warm_kinetic_energy");
  add_entry("ISM_cold_kinetic_energy");

}


void HSTManager::UserWorkBeforeHstOutput(MeshBlock* pmb) {

  auto hst_data = hst->get_proxy(pmb);  // 用于输出到 hst

  // 简写 pmb 的各个成员变量
  Coordinates *pcoord = pmb->pcoord;
  Hydro *phydro = pmb->phydro;
  EquationOfState *peos = pmb->peos;
  Mesh *pmy_mesh = pmb->pmy_mesh;
  const int ks = pmb->ks, ke = pmb->ke, js = pmb->js, je = pmb->je, is = pmb->is, ie = pmb->ie;

  //* 清零：对于「当前时刻的统计量」，应当先清零，再在遍历 kji 时累加。
  // 清零角动量
  for (auto dir : std::string("xyz")) {  // 必须用 std::string("xyz")，否则 "xyz" 是 const char[4]，隐含一个尾随的 '\0'!
    hst_data[std::string("total_angular_momentum_") + dir] = 0.0;
  }
  // 清零 multi-phase ISM 统计量
  for (auto phase : {"hot", "warm", "cold"}) {
    for (auto prop : {"volume", "mass", "thermal_energy", "kinetic_energy"}) {
      hst_data[std::string("ISM_") + phase + "_" + prop] = 0.0;
    }
  }

  AthenaArray<Real> vol(pmb->ncells1);
  for (int k=ks; k<=ke; k++) {
    Real z = pcoord->x3v(k);
    for (int j=js; j<=je; j++) {
      Real y = pcoord->x2v(j);
      pcoord->CellVolume(k, j, is, ie, vol);
      // 这里似乎不能 SIMD，因为 reduction 比较复杂
      for (int i=is; i<=ie; i++) {
        Real x = pcoord->x1v(i);
        Real cell_volume = vol(i);

        // 动量 (用于计算角动量)
        Vector momentum;
        for (int l = 0; l < 3; ++l) {
          momentum[l] = phydro->u(IM1+l,k,j,i) * cell_volume;
        }
        // 总角动量
        Vector angular_momentum = CrossProduct({x,y,z}, momentum);   //* 目前假设 BH 位于 {0,0,0}，否则 {x,y,z} 要减去 BH 的位置，而且速度也要改为相对速度！
        for (int l = 0; l < 3; ++l) {
          hst_data[std::string("total_angular_momentum_") + "xyz"[l] ] += angular_momentum[l];
        }

        // multi-phase ISM 的统计量
        const Real rho = phydro->w(IDN,k,j,i);
        const Real P   = phydro->w(IPR,k,j,i);
        const Real E_thermal = P / (peos->GetGamma() - 1.0);
        const Real E_k       = phydro->u(IEN,k,j,i) - E_thermal;
        const Real T_cgs = Abundance::mu * P / rho * pmy_mesh->punit->code_temperature_mu_cgs;
        
        std::string phase;
        if (T_hot_warm <= T_cgs) { // hot
          phase = "hot";
        } else if (T_warm_cold <= T_cgs and T_cgs < T_hot_warm ) { // warm
          phase = "warm";
        } else { // cold  //BUG 如果是 nan 怎么算？目前这样算成 cold 了。但是如果 phase 是三者之外的某个字符串，那么后续 hst_data 索引会报错！可能用 break？
          phase = "cold";
        }

        hst_data["ISM_" + phase + "_volume"] += cell_volume;
        hst_data["ISM_" + phase + "_mass"] += rho * cell_volume;
        hst_data["ISM_" + phase + "_thermal_energy"] += E_thermal * cell_volume;
        hst_data["ISM_" + phase + "_kinetic_energy"] += E_k * cell_volume;

        //TODO 分不同的 R 内进行统计。目前只是草稿。
        // Real r = std::sqrt(x*x + y*y + z*z);
        // for (Real R: {1,2,3}) {
        //   if (r < R) {
        //     hst_data["ISM_" + phase + "_volume_r<" + std::to_string(R)] += cell_volume;
        //     hst_data["ISM_" + phase + "_mass_r<" + std::to_string(R)] += rho * cell_volume;
        //     hst_data["ISM_" + phase + "_thermal_energy_r<" + std::to_string(R)] += E_thermal * cell_volume;
        //     hst_data["ISM_" + phase + "_kinetic_energy_r<" + std::to_string(R)] += E_k * cell_volume;
        //   }
        // }
      }
    }
  }

  // hst 中的瞬时量
  hst_data["num_MeshBlocks"] = 1.0;
  hst_data["dt_hyperbolic"] = pmy_mesh->dt_hyperbolic;
  hst_data["dt_user"] = pmy_mesh->dt_user;

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


  // 从 input file 中读取参数，存到 pgen 文件的模块级别变量中
  Real M_BH_in_Msun = pin->GetOrAddReal("problem","M_BH", 0.0);
  M_BH = M_BH_in_Msun * punit->solar_mass_code;
  GM_BH = M_BH * punit->grav_const_code;

  R_in = pin->GetOrAddReal("problem","R_in", 0.0);
  R_out = pin->GetOrAddReal("problem","R_out", std::numeric_limits<Real>::max());

  // 读取 BH sink region 的参数，设定 sink 内的密度
  Real T_sink = pin->GetOrAddReal("problem","T_sink", 1.0);
  Real float_min = std::numeric_limits<float>::min();
  Real pressure_floor = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min));
  Real mu = Abundance::mu;
  rho_sink = mu * pressure_floor * punit->hydrogen_mass_code / (punit->k_boltzmann_code * T_sink);

  // hst 中区分不同 phase 的温度，in cgs
  T_hot_warm = pin->GetOrAddReal("hst", "T_hot_warm", 1e5);
  T_warm_cold = pin->GetOrAddReal("hst", "T_warm_cold", 1e3);


  // 构造我自定义机制的对象
  //* 注意：在这里构造的对象，都必须是确定性的，从而保证在 restart 时 / 跨 MPI rank 的一致性。
  initial_condition = std::unique_ptr<InitialCondition>(new InitialCondition(this, pin));

  cooling = Cooling(this, pin);  // 初始化 cooling

  supernovae = Supernovae(this, pin);  // 初始化 supernovae

  progress_report = std::unique_ptr<ProgressReport>(new ProgressReport(this, pin));  // 初始化 progress_report

  hst = std::unique_ptr<HSTManager>(new HSTManager(this, RealUserMeshBlockDataIndex::hst));  // 初始化 hst 机制

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

  // 手动给各个 Passive Scalar 场起名，输出到 info 下面的文件，从而指导后处理读取相应的 output （包括 hst 中的 %d-scalar）
  // 目前先 hard code，每当有更改时要手动更新。未来可以考虑从 input file 中读取。
  std::map<PassiveScalarIndex::PassiveScalarIndex, std::string> passive_scalar_names;
  passive_scalar_names[PassiveScalarIndex::SN]             = "SN_tracer";
  passive_scalar_names[PassiveScalarIndex::initial_radius] = "initial_radius";
  passive_scalar_names[PassiveScalarIndex::initial_fluid ] = "initial_fluid";


  if (Globals::my_rank == 0) {
    std::cout << "pgen 编译于: " << __DATE__ << " " << __TIME__ << '\n'; 
    // Print unit 相关信息
    punit->PrintCodeUnits();
    punit->PrintConstantsInCodeUnits();
    // 把 passive scalar 的名字写入文件，从而指导后处理在 output 文件中读取相应的 field
    writeMapToFile(passive_scalar_names, "info/passive_scalar.txt");
    //TODO Print 自定义的参数
  }

  return;
}

// 调用时机：程序启动时（包括 restart 时）
// 函数用途：初始化每个 MeshBlock 上的局部变量（包括数组数据）
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // 初始化 ruser_meshblock_data
  { // 防止 idx 别名泄露
    namespace idx = RealUserMeshBlockDataIndex;
    AllocateRealUserMeshBlockDataField(idx::N_RealUserMeshBlockData);
    hst->RequestMeshBlockData(this);
  }

  // 设定 UOV
  {
    using namespace UOV;
    AllocateUserOutputVariables(N_UOV);  // allocate 我自定义的 output 变量（uov, user_out_var）
    //TODO 把 enum 改为使用 vector 或 map 之类的，从而根据开启什么功能来决定总共 UOV 的数目，减小输出体积。但是，从代称到编号的 map 需要在其他文件中也可见！
    //* 为每个 user_out_var 变量设置名字。这里使用 _snake_case，额外在开头加一个 _ 以规避 yt 中名称重复导致的不便。
    if (cooling.cooling_flag) {
      SetUserOutputVariableName(cooling_rate, "_cooling_rate");
    }
  }


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

  if (progress_report->is_time_to_report()) {
    progress_report->report();
  }
  return;
}


// 调用时机：每个将要输出 output 的时间步的末尾。(但不包括 hst 输出，因为太频繁了) 
//* 生在 UserWorkInLoop 之后，甚至在 pmesh->time 和 pmesh->dt 已经更新后
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

  // 在结束时，最后一次 progress_report
  progress_report->report(true);
}
