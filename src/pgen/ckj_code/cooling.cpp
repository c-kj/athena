
// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"
#include "../../mesh/mesh.hpp"

// 自定义的头文件
#include "cooling.hpp"

void CoolingSourceTerm(std::string cooling_model,  //* 这是临时处理，把 cooling_type 作为参数传递进来
            MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar) {

  auto punit = pmb->pmy_mesh->punit;
  Real mu = 1.0;  //? mu 取什么值？
  Real gamma = pmb->peos->GetGamma();

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real cooling_func_value_cgs = 0.0; // 必须初始化，否则后面 if 分支如果没有赋值的话，是 UB！会随机给一个值。

        Real rho = prim(IDN,k,j,i); 
        Real P = prim(IPR,k,j,i);
        Real T_cgs = mu * P / rho * punit->code_temperature_mu_cgs;
        Real n_cgs = rho * punit->code_density_cgs / (mu * Constants::hydrogen_mass_cgs); 
        Real dt_cgs = dt * punit->code_time_cgs;
        
        if (cooling_model == "supernova") { 
          // Drain 2011, ISM textbook, Sec 34.1, eqn 34.2 & 34.3。用于估计 SN shock cooling 的简单幂律函数。高于 10^(7.3) K 的部分是轫致辐射主导的。
          if (/*  1e5 < T_cgs &&  */ T_cgs <= std::pow(10.0, 7.3) ){
            cooling_func_value_cgs = 1.1e-22 * std::pow(T_cgs/1e6, -0.7) * std::exp(-1.18348e5/T_cgs); //TEMP 指数截断
          } else if ( std::pow(10.0, 7.3) < T_cgs ) {
            cooling_func_value_cgs = 2.3e-24 * std::pow(T_cgs/1e6, 0.5);
          }
          Real dE = - cooling_func_value_cgs * n_cgs * n_cgs * dt_cgs / punit->code_energydensity_cgs;
          dE = std::max(dE, -P/(gamma-1) * 0.1);  //TEMP 简陋的限制：让一个 timestep 内热能的减少量不会超过原来热能的 10%。比例越大，限制越弱。
          // Real new_T_cgs = mu * (P - dE*(gamma-1)) / rho * punit->code_temperature_mu_cgs;
          Real original_energy_density = cons(IEN,k,j,i); // 只在 debug 时用来查看原来的 E 的值

          cons(IEN,k,j,i) += dE;
        } else {
          throw std::invalid_argument("cooling_type not supported");
        }
      }
    }
  }
  return;
}
