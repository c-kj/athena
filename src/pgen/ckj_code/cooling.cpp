
// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"

// 自定义的头文件
#include "cooling.hpp"


Cooling::Cooling(ParameterInput *pin, Mesh *pmy_mesh): punit(pmy_mesh->punit) {  // 在初始化列表中把传入的 punit 赋值给成员变量 punit

  cooling_model = pin->GetOrAddString("cooling", "cooling_model", "none");
  cooling_flag = cooling_model != "none";  // 如果没有指定 cooling_model 或压根没有 cooling 这个 block，则 cooling_flag 为 false
  if (!cooling_flag) return;  // 如果没有开启 cooling，则不继续初始化
  
  CFL_cooling = pin->GetOrAddReal("cooling", "CFL_cooling", 1.0); 
  if (CFL_cooling <= 0.0) {throw std::invalid_argument("CFL_cooling must be positive!");}
  // operator_splitting = pin->GetOrAddBoolean("cooling", "operator_splitting", true);
  implicit = pin->GetOrAddBoolean("cooling", "implicit", true);
  integrator = pin->GetOrAddString("cooling", "integrator", "Euler");


  // 设定元素丰度 //* 目前暂时放在 Cooling 类的成员变量中，但如果组分要演化，则再考虑更改。
  X_H = 0.7, X_Metal = 0.01295;  // H, He, Metal 元素的质量分数
  // X_H = 1.0, X_Metal = 0.0;   // 纯 H
  X_He = 1.0 - X_H - X_Metal;

  mu = 4.0/(5*X_H + 3 - X_Metal);  //? mu 取什么值？
  mu_e = 2.0/(1+X_H);


}


// 返回冷却函数的值 (in cgs unit)
const Real Cooling::CoolingFunction(Real T_cgs) {
  if (T_cgs <= 0.0) return 0.0; // 避免 T<0 时指数截断反而变得巨大
  if (cooling_model == "supernova") { 
  // Drain 2011, ISM textbook, Sec 34.1, eqn 34.2 & 34.3。用于估计 SN shock cooling 的简单幂律函数。高于 10^(7.3) K 的部分是轫致辐射主导的。
    if (/*  1e5 < T_cgs &&  */ T_cgs <= std::pow(10.0, 7.3) ){
      return 1.1e-22 * std::pow(T_cgs/1e6, -0.7) * std::exp(-1.18348e5/T_cgs); //TEMP 指数截断，具体的值有待调整
    } else if ( std::pow(10.0, 7.3) < T_cgs ) {
      return 2.3e-24 * std::pow(T_cgs/1e6, 0.5);
    }
  } else {
    throw std::invalid_argument("cooling_model not supported");
  }
  return 0.0;  // 默认返回值
}


// 计算 Cooling Rate (in code unit)，即 d(能量密度) / dt。//* 注意：这里 Cooling Rate 是正数！
const Real Cooling::CoolingRate(const Real &rho, const Real &P) { //TODO 这里有问题：使用 RK4 的时候会需要多次传入 T ？
  Real T_cgs = mu * P / rho * punit->code_temperature_mu_cgs;

  // 计算数密度 n
  Real rho_cgs = rho * punit->code_density_cgs;
  // Real n_cgs = rho_cgs / (mu * Constants::hydrogen_mass_cgs); 
  Real n_e_cgs = rho_cgs / (mu_e * Constants::hydrogen_mass_cgs);
  Real n_H_cgs = rho_cgs * X_H / Constants::hydrogen_mass_cgs;  //? 这么算适用于 Draine 2011 吗？ 适用于其他 Cooling Function 吗？
  
  //? 使用 n_e * n_H 还是 n * n ? 这取决于文献给出 Cooling Function 时如何做归一化。
  Real cooling_rate_code = n_e_cgs * n_H_cgs * CoolingFunction(T_cgs) / (punit->code_energydensity_cgs /  punit->code_time_cgs); 
  
  return cooling_rate_code;
}


const Real Cooling::CoolingTimeScale(const Real &E_thermal, const Real &cooling_rate) {
  if (cooling_rate == 0.0) return std::numeric_limits<Real>::max();  // 当 cooling_rate 为 0 时，冷却时标设为一个很大的数。
  //? 要处理 E_thermal <= 0 的情况吗？
  return E_thermal / cooling_rate;
}


const Real Cooling::CoolingTimeStep(MeshBlock *pmb) {
  Real cooling_dt = std::numeric_limits<Real>::max();  // 先初始化为一个很大的数
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {  // 这里涉及归约操作，Athena++ 的 OMP 好像是基于 MeshBlock 的，搞不清楚，这里就不用 OMP 了
        Real rho = pmb->phydro->w(IDN,k,j,i);
        Real P = pmb->phydro->w(IPR,k,j,i); 

        //? 这里先计算 dt_cooling ，后续又调用 CoolingRate，重复计算了。可以考虑优化？但如果是 RK4 之类的积分器，其实总共要调用 4 次（各不同），可能意义不大？
        Real cooling_rate = CoolingRate(rho, P);
        Real E_thermal = P / (pmb->peos->GetGamma() - 1.0);
        Real dt = CFL_cooling * CoolingTimeScale(E_thermal, cooling_rate);   //calculate your own time step here
        cooling_dt = std::min(cooling_dt, dt);
      }
    }
  }
  return cooling_dt;
}


void Cooling::CoolingSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar) {

  Real T_floor_cgs = 1e2; // TEMP

  Real gamma = pmb->peos->GetGamma(); //TEMP 只适用于理想气体 EoS 的情况


  // Real t_cool = std::numeric_limits<Real>::max(); 

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {  // Athena++ 的 OMP 好像是基于 MeshBlock 的，搞不清楚，如果涉及归约操作的话，OMP 并行可能有问题。先不启用。

        Real rho = prim(IDN,k,j,i); 
        Real P = prim(IPR,k,j,i);
        //TEMP 尝试使用 cons 而非 prim？好像并不合理
        // Real rho = cons(IDN,k,j,i); 
        // Real P = (gamma - 1.0) * (cons(IEN,k,j,i) - 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i))) / rho);  // 用 EoS 计算 P
        if (P <= 0.0) { //TEMP
          std::cout << "P = " << P << std::endl; 
          throw std::invalid_argument("P <= 0.0");
        }
        
        // Real E_thermal = pmb->peos->EgasFromRhoP(rho, P);  //BUG 现在不知道为啥有 bug
        Real E_thermal = P / (gamma - 1.0);  //* 只适用于 Ideal Gas
        //BUG 这里似乎不应该用 prim 来算，而是从 cons 中算出 E_thermal，这样才是 Limiter 所需要的「当前剩下的 E_thermal」

        //TEMP Integrator，待整理
        Real dE = 0.0;
        if (!implicit) {  // 显式 integrator
          if (integrator == "Euler") {  // Forward Euler
            Real cooling_rate = CoolingRate(rho, P);
            dE = - cooling_rate * dt; 
          } else {
            throw std::invalid_argument("integrator not supported");
          }
        } else {  // 隐式 integrator
          if (integrator == "Euler") {  // Implicit Euler 算法, 不动点迭代
            Real cooling_rate = 0.0;
            for (int l = 0; l < 5; l++) {
              cooling_rate = CoolingRate(rho, P + dE*(gamma - 1.0));
              dE = - cooling_rate * dt;
            }
          } else if (integrator == "Semi-Implicit Euler") {  // Semi-Implicit Backward Euler 算法，这个命名是暂时的
            Real dt_subcycle = std::min(CoolingTimeScale(E_thermal, CoolingRate(rho, P)) * CFL_cooling / 5.0, dt); //TEMP
            Real dE_subcycle = 0.0;
            for (Real t=0.0; t<dt; t+=dt_subcycle) {
              P += dE_subcycle*(gamma - 1.0);
              Real cooling_rate = CoolingRate(rho, P);
              Real dP = 1e-3 * P;  //TEMP 用有限差分计算 Jacobian
              Real jacobian = (gamma-1.0)*((-CoolingRate(rho,P+dP)) - (-cooling_rate)) / (dP); // 注意 RHS 是 - CoolingRate
              dE_subcycle = dt_subcycle / (1 - jacobian*dt_subcycle) * (-cooling_rate);
              dE += dE_subcycle;
            }
          } else {
            throw std::invalid_argument("integrator not supported");
          }
        }

        // 只用来诊断
        // t_cool = std::min(CoolingTimeScale(E_thermal, cooling_rate), t_cool); //TEMP

        // dE = - std::min(std::abs(dE), E_thermal * 0.1);  //TEMP 简陋的限制：让一个 timestep 内热能的减少量不会超过原来热能的 10%。比例越大，限制越弱。
        // Real new_T_cgs = mu * (P - dE*(gamma-1)) / rho * punit->code_temperature_mu_cgs;

        //TEMP 限制最低温度。但要考虑到 source term 在一个 dt 内被调用多次
        Real E_thermal_floor = rho * T_floor_cgs / mu / punit->code_temperature_mu_cgs / (gamma - 1);
        dE = - std::min(std::abs(dE),  (E_thermal-E_thermal_floor));  


        cons(IEN,k,j,i) += dE;
      }
    }
  }
  return;
}