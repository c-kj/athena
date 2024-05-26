#pragma once

// 用于实现初始条件。同时提供一个 SetSingleCell 函数，用于在源项中将「外部区域」（r > R_out）的流体状态固定为初始条件不变。
// 目前的实现思路是：在 InitUserMeshData 中初始化这个类的实例，读入必要的参数；然后在 SetSingleCell 的具体分支中，读入各种 init_cond_type 所需要的参数，用 static 来使得只读取一次。

// C++ headers
#include <string>  // std::string

// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../hydro/hydro.hpp"

// 自定义的头文件
#include "ckj_code.hpp"
#include "region.hpp"

#if GENERAL_EOS  // 目前不支持 General EOS，因为从输入的 n、T 转换到 code 所需的 rho、E_thermal 需要更多考量
  #error "General EOS is not supported yet for init_condition.hpp !"
#endif

struct InitialCondition {
  // 需要的指针
  ParameterInput *pin;
  Units *punit;
  // 读入的参数
  const Real gamma;
  const std::string init_cond_type;
  // 从 T 转换到 P 需要平均分子量 mu，这里暂时 hard code 为常量
  static constexpr Real mu = Abundance::mu; 
  // 初始条件中最根本的参数（总是需要）
  const Real n_init_cgs;
  const Real T_init_cgs;
  // 换算出 code unit 下的 rho 和 E_thermal
  Real rho_init_code;
  Real E_thermal_init_code;

  InitialCondition(Mesh *pmesh, ParameterInput *pin);
  // 成员函数
  void SetInitialCondition(MeshBlock *pmb);
  void SetSingleCell(MeshBlock *pmb, const int i, const int j, const int k);
  void PrintInfo();  //TODO 在初始化时打印初值的信息

  // 用于 approximate_Bondi 的 inline 函数
  Real approx_Bondi_rho_profile(Real alpha, Real R_Bondi, Real r) {
    return pow(1 + alpha * R_Bondi / r, 1.5);
  }

};

// 在 cgs unit 下，将 n 换算为 rho
inline Real rho_from_n_cgs(Real mu, Real n_cgs) {
  Real rho_cgs = mu * n_cgs * Constants::hydrogen_mass_cgs;
  return rho_cgs;
}

// 在 cgs unit 下，将 n 和 T 换算为 P
inline Real P_from_nT_cgs(Real n_cgs, Real T_cgs) {
  Real P_cgs = n_cgs * Constants::k_boltzmann_cgs * T_cgs;
  return P_cgs;
}
