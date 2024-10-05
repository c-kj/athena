#pragma once

#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"


#include "my_outputs.hpp"

//* 添加新 hst output 的步骤
// 1. 在 hst_index::UserHistoryOutputIndex 中添加新的 index
// 2. 在 hst_funcs::hst_XXX 中添加新的函数
// 3. 在 MeshBlock::UserWorkBeforeOutput 中 EnrollUserHistoryOutput
// 如果需要累积储存，还需要配合 MeshBlock 上的数据：
// 1. 在RealUserMeshBlockDataIndex 中添加新的 index
// 2. MeshBlock::InitUserMeshBlockData 中初始化 ruser_meshblock_data
// 3. 在相应的源项中修改 ruser_meshblock_data 的值


// 自定义的 hst Output 的 index。用于 EnrollUserHistoryOutput 时；可能用于 hst_ 系列函数中，与 iout 参数做判断。
namespace hst_index {
  enum UserHistoryOutputIndex {
    num_MeshBlocks,
    dt_hyperbolic,
    dt_user,
    // 源项造成的改变量
    // 吸积量
    accreted_mass,
    accretion_rate,
    accreted_SN_tracer,
    accretion_rate_SN_tracer,
    accreted_energy,
    accreted_momentum_x,
    accreted_momentum_y,
    accreted_momentum_z,
    accreted_angular_momentum_x,
    accreted_angular_momentum_y,
    accreted_angular_momentum_z,
    // 能量
    total_cooling_loss,
    BH_gravity_work,
    //TODO BH_potential_energy,  等把 BH 的相关信息抽出来作为全局变量时再写
    // SN 注入量
    SN_injected_energy,
    SN_injected_mass,
    SN_injected_number,
    // 其他 ...
    N_UserHistoryOutput // 总个数
  };
}


/* -------------------------------------------------------------------------- */
/*                          自定义 hst Output 函数                              */
/* -------------------------------------------------------------------------- */

namespace hst_funcs {

namespace idx = RealUserMeshBlockDataIndex; // 这里只用得到 RealUserMeshBlockData，因此用别名来简化名称

#if false  // 发现只要不 Enroll，就会输出 0，并且 header 的名字为空。因此目前不需要这个函数
// hst_NaN 用于在某些配置下，直接输出 nan。这样在读取时就可以很容易发现相应的 hst 输出是不启用的。
// 在 EnrollUserHistoryOutput 时，使用三目运算符来使用 hst_NaN，这样应该比在相应的 hst function 里返回 NaN 要更快一点。
Real hst_NaN(MeshBlock *pmb, int iout) {
  return std::numeric_limits<Real>::quiet_NaN();
}
#endif

// 统计总共有多少个 MeshBlocks
Real hst_num_MeshBlocks(MeshBlock *pmb, int iout) {
  return 1;  // 每个 MeshBlock 返回 1，最后加起来。
}

Real hst_dt_hyperbolic(MeshBlock *pmb, int iout) {
  return pmb->pmy_mesh->dt_hyperbolic;
}

Real hst_dt_user(MeshBlock *pmb, int iout) {
  return pmb->pmy_mesh->dt_user;  
  // 如果没 EnrollUserTimeStepFunction，这里的 dt_user 会是 Athena++ 内部初始化时设定的 std::numeric_limits<Real>::max()
}

// 累计吸积的质量
Real hst_accreted_mass(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_mass](0);
}
// 输出当前时间步内的 accretion rate
Real hst_accretion_rate(MeshBlock *pmb, int iout) {
  Real accretion_rate = pmb->ruser_meshblock_data[idx::accretion_rate](0);
  return accretion_rate;
}

// 累计吸积的 SN tracer 质量
Real hst_accreted_SN_tracer(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_SN_tracer](0);
}
// 输出当前时间步内的 accretion rate of SN tracer
Real hst_accretion_rate_SN_tracer(MeshBlock *pmb, int iout) {
  Real accretion_rate_SN_tracer = pmb->ruser_meshblock_data[idx::accretion_rate_SN_tracer](0);
  return accretion_rate_SN_tracer;
}

// 累计吸积的能量
Real hst_accreted_energy(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_energy](0);
}

// 动量的三个分量
Real hst_accreted_momentum_x(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_momentum](0);
}
Real hst_accreted_momentum_y(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_momentum](1);
}
Real hst_accreted_momentum_z(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_momentum](2);
}

// 角动量的三个分量
Real hst_accreted_angular_momentum_x(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_angular_momentum](0);
}
Real hst_accreted_angular_momentum_y(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_angular_momentum](1);
}
Real hst_accreted_angular_momentum_z(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::accreted_angular_momentum](2);
}

// 源项的改变量
Real hst_total_cooling_loss(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::total_cooling_loss](0);
}

Real hst_BH_gravity_work(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::BH_gravity_work](0);
}

// Real hst_BH_potential_energy(MeshBlock *pmb, int iout) {
//   for (int k = pmb->ks; k <= pmb->ke; ++k) {
//     Real z = pmb->pcoord->x3v(k);
//     for (int j = pmb->js; j <= pmb->je; ++j) {
//       Real y = pmb->pcoord->x2v(j);
// #pragma omp simd  // SIMD 矢量化。不确定是否成功。
//       for (int i = pmb->is; i <= pmb->ie; ++i) {
//         Real x = pmb->pcoord->x1v(i);
//         Real r = std::sqrt(x*x + y*y + z*z);
//         Real r3 = r*r*r;
//         Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);

//         Real dens = pmb->phydro->u(IDN,k,j,i); // 用 dens 称呼 cons 守恒密度，与 prim 密度区分。

//       }
//     }
//   }
// }

Real hst_SN_injected_energy(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::SN_injected_energy](0);
}
Real hst_SN_injected_mass(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::SN_injected_mass](0);
}
Real hst_SN_injected_number(MeshBlock *pmb, int iout) {
  return pmb->ruser_meshblock_data[idx::SN_injected_number](0);
}


} // namespace hst_funcs