
// Athena++ headers
#include "../../eos/eos.hpp"
#include "../../hydro/hydro.hpp"

// 自定义的头文件
#include "cooling.hpp"
#include "my_outputs.hpp"
#include "ckj_plugin.hpp"


Cooling::Cooling(Mesh *pmy_mesh, ParameterInput *pin): punit(pmy_mesh->punit) {  // 在初始化列表中把传入的 punit 赋值给成员变量 punit

  //TODO 要不要用初始化列表？
  // 使用初始化列表的主要好处：可以初始化 const 和 引用 类型的成员变量

  // 配置 CoolingModel & 是否开启 cooling
  cooling_model = pin->GetOrAddString("cooling", "cooling_model", "none");
  cooling_flag = cooling_model != "none";  // 如果没有指定 cooling_model 或压根没有 cooling 这个 block，则 cooling_flag 为 false
  if (!cooling_flag) {return;}  // 如果没有开启 cooling，则不继续初始化
  model = CoolingModel::Create(cooling_model);  // 根据输入的 cooling_model 字符串，创建一个 CoolingModel 的指针
  
  // 从 input 中确定 CoolingSourceTerm 在什么时机注入
  // 放在 SourceTerm 与 UserWorkInLoop 的区别是：SourceTerm 中 Cooling 使用的是未被其它 SourceTerm（比如 SN、Grav）修改的 prim。（若使用 cons 则是被修改过的）
  // 如果这一步注入 SN，那么 SourceTerm 并不会在这一步就产生强冷却，但 UserWorkInLoop 会立刻产生强冷却（而且步长未被限制）。
  //* 目前看来似乎 InSourceTerm 比较好（虽然要算两次）。如果后续确定下来，可以改为在 .hpp 中以 static constexpr 定义，从而在编译器优化掉。
  std::string source_term_position_str = pin->GetOrAddString("cooling", "source_term_position", "InSourceTerm");
  source_term_position = source_term_position_map.at(source_term_position_str);
  use_prim_in_cooling = pin->GetOrAddBoolean("cooling", "use_prim_in_cooling", true);  // 使用 prim 或 cons 来计算 rho、P


  // 配置 integrator
  integrator = pin->GetOrAddString("cooling", "integrator", "Euler");
  implicit = pin->GetOrAddBoolean("cooling", "implicit", true);
  CFL_cooling = pin->GetOrAddReal("cooling", "CFL_cooling", 1.0); 
  if (CFL_cooling <= 0.0) {throw std::invalid_argument("CFL_cooling must be positive!");}

  // 配置 RootFinder
  max_iter = pin->GetOrAddInteger("cooling", "max_iter", 10);  // 默认 10 次迭代，应该足够了。有 rel_tol 存在，一般应该不至于需要 10 次。
  rel_tol = pin->GetOrAddReal("cooling", "rel_tol", 1e-6);     //? 多少合适？
  abs_tol = pin->GetOrAddReal("cooling", "abs_tol", 0);        // 默认绝对容差为 0，因为这个值不太好把握，容易造成一次迭代都没有就直接退出。

  // 配置 Subcycle
  CFL_subcycle = pin->GetOrAddReal("cooling", "CFL_subcycle", 1.0);
  subcycle_adaptive = pin->GetOrAddBoolean("cooling", "subcycle_adaptive", true);

  // 配置 Limiter
  limiter_on = pin->GetOrAddBoolean("cooling", "limiter_on", false);  // cooling 的结尾是否要用 limiter 限制
  T_floor_cgs = pin->GetOrAddReal("cooling", "T_floor", 1e2);

}




// 计算 Cooling Rate (in code unit)，即 d(能量密度) / dt。//* 注意：这里 Cooling Rate 是正数！
Real Cooling::CoolingRate(const Real rho, const Real P) const {
  Real T_cgs = mu * P / rho * punit->code_temperature_mu_cgs;

  // 计算数密度 n
  Real rho_cgs = rho * punit->code_density_cgs;
  //TODO 检查这里的公式对吗？
  // Real n_cgs = rho_cgs / (mu * Constants::hydrogen_mass_cgs); 
  Real n_e_cgs = rho_cgs / (mu_e * Constants::hydrogen_mass_cgs);
  Real n_H_cgs = rho_cgs * X_H / Constants::hydrogen_mass_cgs;  //? 这么算适用于 Draine 2011 吗？ 适用于其他 Cooling Function 吗？
  
  //? 使用 n_e * n_H 还是 n * n ? 这取决于文献给出 Cooling Function 时如何做归一化。
  Real cooling_rate_code = n_e_cgs * n_H_cgs * model->CoolingFunction(T_cgs) / (punit->code_energydensity_cgs /  punit->code_time_cgs); 
  
  return cooling_rate_code;
}


// 计算 Cooling 时标（不乘额外系数）。返回值一定大于 0 （如果 < 0 则返回正无穷大，从而不对 TimeStep 做任何限制）
Real Cooling::CoolingTimeScale(const Real E_thermal, const Real cooling_rate) {
  if (cooling_rate == 0.0) { return std::numeric_limits<Real>::max(); }  // 当 cooling_rate 为 0 时，冷却时标设为一个很大的数，避免除 0 错误
  Real cooling_timescale = E_thermal / cooling_rate;
  if (cooling_timescale <= 0.0) { return std::numeric_limits<Real>::max(); }  // 如果得到的 cooling_timescale <= 0，则返回一个很大的数，从而保证返回值总是 > 0 的
  return cooling_timescale;
}

Real Cooling::CoolingTimeStep(MeshBlock *pmb) const {
  Real gm1 = (pmb->peos->GetGamma() - 1.0);
  Real cooling_dt = std::numeric_limits<Real>::max();  // 先初始化为一个很大的数
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
// #pragma omp simd reduction(min:cooling_dt)  // SIMD 矢量化加速。似乎因为内部调用太复杂无法内联，目前无法矢量化。目前不知道这里是不是性能热点，暂时不用。
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rho = pmb->phydro->w(IDN,k,j,i);
        Real P = pmb->phydro->w(IPR,k,j,i); 

        //? 这里先计算 dt_cooling ，后续又调用 CoolingRate，重复计算了。可以考虑优化？但如果是 RK4 之类的积分器，其实总共要调用 4 次（各不同），可能意义不大？
        Real cooling_rate = CoolingRate(rho, P);
        Real E_thermal = P / gm1;
        Real dt = CFL_cooling * CoolingTimeScale(E_thermal, cooling_rate);   //calculate your own time step here

        cooling_dt = std::min(cooling_dt, dt);
      }
    }
  }
  return cooling_dt;
}

// 求解 dy/dt = RHS(y)，给出 dy
// 目前这里的 y 指代的是 E_thermal
// 这里假定了 RHS 不是 t 的函数。对于 Cooling 来说，总是成立的。
Real Cooling::Integrator(const std::function<Real(Real)>& RHS, Real y0, Real dt) const {
  Real dy = 0.0;
  if (!implicit) { // 显式 integrator
    if (integrator == "Euler") {
      dy = RHS(y0) * dt;
    } else if (integrator == "RK4") {
      Real k1 = RHS(y0);
      Real k2 = RHS(y0 + 0.5 * k1 * dt);
      Real k3 = RHS(y0 + 0.5 * k2 * dt);
      Real k4 = RHS(y0 + k3 * dt);
      dy = (k1 + 2*k2 + 2*k3 + k4) * dt / 6.0;
    } else {
      throw std::invalid_argument("integrator not supported");
    }
  } else { // 隐式 integrator
    std::function<Real(Real)> func;  // 若为隐式离散化，则把方程整理成关于待求的 dy 的函数 func，用 RootFinder 求根
    if (integrator == "Euler") {
      func = [&RHS,y0,dt](Real dy) {return RHS(y0+dy)*dt - dy;};  // 用于RootFinder 的函数
    } else {
      throw std::invalid_argument("integrator not supported");
    }

    dy = RootFinder(func, 0.0, max_iter, rel_tol, abs_tol); // 这里对 dy 的初始猜测设为 0 。不能由 Forward Euler 给出，因为 Forward Euler 可能上来就过冲了。
  }
  return dy;
}

// 求解 func(x) == 0 的方程，给出 x 的根。初始猜测为 x0
// 目前这里的 x 指代的是 dy，即 d E_thermal
Real Cooling::RootFinder(const std::function<Real(Real)>& func, Real x0, int max_iter, Real rel_tol, Real abs_tol) {
  Real x = x0;
  if (func(x) == 0.0) {return x;}  // 如果 Cooling Function 给出的 dE == 0，那么无需求根，直接返回。这也避免了对于临近 x，func 都为 0 时会导致的 0/0 问题。
  for (int i = 0; i < max_iter; ++i) {  // 最多迭代 max_iter 次。若 max_iter 为 1，则为 ODE 的半隐式方法
    //* 目前，在进行迭代前进行容差判断。坏处：如果容差过大（abs_tol 不好把握），可能连一次半隐式迭代都没有，得到虚假的 dE = 0.
    Real tol = rel_tol * std::abs(x) + abs_tol;  // 在每次迭代中，都根据 x 的值重新计算 tol  //TODO 这里的容差究竟应该如何计算，有待重新考虑。
    if (std::abs(func(x)) <= tol) { break; } // 如果满足精度要求，则退出循环

    // Newton-Raphson 迭代
    Real df_dx = FiniteDifferenceDerivative(func, x);  // 目前使用有限差分来计算导数
    x = x - func(x) / df_dx;  //? 如果 df_dx == 0 怎么办？
  }
  return x;
}

// 求函数 func(x) 在 x 处的导数值
Real Cooling::FiniteDifferenceDerivative(const std::function<Real(Real)>& func, Real x) {
  // 这部分参照 Numerical Recipes 5.7 的内容
  static const Real epsilon = std::numeric_limits<Real>::epsilon();
  Real h = std::pow(epsilon, 1.0/3.0);  // 对于二阶精度中心差分，幂次是 1/3; 如果是一阶精度单侧差分，幂次要改为 1/2
  h *= std::max(1.0, std::abs(x));  // 乘以相应的特征长度 x，但若 x 太小则截断为 1.0

  // 调整 h 的值以尽量避免舍入误差
  volatile Real x_plus_h = x + h;
  h = x_plus_h - x;  // 重新计算 h，以确保 h 是可被浮点数精确表示

  return (func(x + h) - func(x - h)) / (2*h); // 使用中心差分，二阶精度
  // 虽然中心差分需要额外计算一次 func 函数，但对于我这里的应用，CoolingFunction 的计算应该代价并不高，而更准确的导数会让 Newton RootFinder 更准确/更快收敛，所以暂时不需要使用单侧差分。
}



void Cooling::CoolingSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar) {

  const Real gamma = pmb->peos->GetGamma(); // 只适用于 adiabatic EoS 的情况，如果是 GENERAL_EOS，则会直接报错。若要使用 GENERAL_EOS，所有调用 gamma 的地方都需要修改！

  //? 这里是否需要避免黑洞内被 Cooling 影响？（至于 R_out 之外应该不用避免 Cooling，否则压强反而不正确了）
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
// #pragma omp simd // SIMD 矢量化加速。如果涉及相邻 cell 的平均，SIMD 可能会有问题。目前不知道这里是不是性能热点，暂时不用。循环内部太复杂，很可能无法矢量化。
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        Real rho, P;
        if (use_prim_in_cooling) {  // 使用 prim 或 cons 来计算 rho、P
          rho = prim(IDN,k,j,i); 
          P = prim(IPR,k,j,i);
        } else {
          rho = cons(IDN,k,j,i); 
          P = (gamma - 1.0) * (cons(IEN,k,j,i) - 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i))) / rho);  
        }

        auto RHS = [gamma,rho,this](Real E_thermal) {
          Real P = (gamma - 1.0) * E_thermal;
          return - CoolingRate(rho, P);  // 注意这里是负号，因为 CoolingRate 是正数，而 dE/dt = RHS 是负数
        };

        Real E_thermal = P / (gamma - 1.0);
        //BUG 这里似乎不应该用 prim 来算，而是从 cons 中算出 E_thermal，这样才是 Limiter 所需要的「当前剩下的 E_thermal」。对于 CoolingTimescale 的计算怎么办？


        Real dE = 0.0;  // 用于最后返回的 dE （整个 SourceTerm 的 dt 内，当前 cell 中的热能密度变化量）。//* 注意：E 是能量密度，而非能量！而且 dE < 0 代表 cooling ！
        Real dt_subcycle = dt;  // 这里初始化的值实际上无所谓，因为在子循环的开始一定会被重新赋值
        Real t_subcycle = 0.0;  // 从 0 增加到 dt
        Real dE_subcycle = 0.0;  // 在当前子循环中增加的 dE
        while (t_subcycle < dt) { // 进入子循环
          if (subcycle_adaptive || t_subcycle == 0.0) {  // 如果 subcycle_adaptive 则每个 subcycle 都重新计算 dt_subcycle，否则只在第一个子循环计算一次
            dt_subcycle = CoolingTimeScale(E_thermal, CoolingRate(rho, P)) * CFL_subcycle; 
          }
          if (t_subcycle + dt_subcycle > dt) {dt_subcycle = dt - t_subcycle;}  // 限制 dt_subcycle 使得 t_subcycle 不能超过 dt，使最后一个子循环结束时 t_subcycle == dt。同时也限制了 dt_subcycle <= dt
          dE_subcycle = Integrator(RHS, E_thermal+dE, dt_subcycle);
          dE += dE_subcycle;
          t_subcycle += dt_subcycle;
        }

        if (limiter_on) {
          // dE = - std::min(std::abs(dE), E_thermal * 0.1);  //TEMP 简陋的限制：让一个 timestep 内热能的减少量不会超过原来热能的 10%。比例越大，限制越弱。
          
          // 限制最低温度
          Real E_thermal_floor = rho * T_floor_cgs / mu / punit->code_temperature_mu_cgs / (gamma - 1.0);
          // 避免 cooling 到 floor 以下，如果本来就在 floor 以下则不管（膨胀冷却、高马赫数伪冷却）
          Real E_thermal_cons = use_prim_in_cooling ? (cons(IEN,k,j,i) - 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i))) / rho) : E_thermal; // 如果前面已经用了 cons 来计算 E_thermal，就不用再算一遍了
          Real E_thermal_minus_floor = E_thermal_cons - E_thermal_floor;    // 当前 E_thermal 与 floor 的差值。不论前面的计算用 prim 还是 cons，这里都用 cons，因为这才是要被改变的量。
          if (E_thermal_minus_floor > 0 and dE < - E_thermal_minus_floor) { // 如果 E_thermal 高于 floor 且 cooling 会把 E_thermal 降到 floor 以下
            dE = - E_thermal_minus_floor;                                   // 则把 dE 设为 - E_thermal_minus_floor，使得 E_thermal 降到 floor 处
          }
        }


        cons(IEN,k,j,i) += dE;
        pmb->user_out_var(UOV::cooling_rate, k,j,i) = -dE/dt;  // 把每个 cell 内的 cooling_rate 储存到 user_out_var 中。
        //* 由于每个 cycle 内，后面的 stage 这里的赋值会覆盖前面的，所以最后只会保留最后一个 stage 的值。不过既然是 rate，那么用最后一个 stage 的值也基本上是准确的。
        // 这里储存 dE/dt 而非 dE 是因为 dE/dt 才是物理的，不涉及人为的 dt。而且如果放在 InSourceTerm 中，那么会被调用半步。
        // 目前储存的是 -dE/dt。若出现加热（由于 Heating Term 或者由于 floor），则可能需要重新考虑应该储存什么

        const Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);
        namespace idx = RealUserMeshBlockDataIndex;
        pmb->ruser_meshblock_data[idx::total_cooling_loss](0) += -dE * cell_volume * ckj_plugin::source_term_weight;  // 历史累计的 cooling loss
      }
    }
  }
  return;
}

// 工厂函数，根据输入的字符串返回对应的 CoolingModel 指针
std::unique_ptr<CoolingModel> CoolingModel::Create(const std::string &cooling_model)
{
  if (cooling_model == "Draine_2011") {
    return std::unique_ptr<CoolingModel>(new Draine_2011());
  } else if (cooling_model == "Draine_2011_cutoff") {
    return std::unique_ptr<CoolingModel>(new Draine_2011_cutoff());
  } else if (cooling_model == "other_model") {
    return std::unique_ptr<CoolingModel>(new OtherModel());
  }
  // 更多模型...
  else {
    throw std::invalid_argument("未定义的 CoolingModel: " + cooling_model);
  }
}

// Drain 2011, ISM textbook, Sec 34.1, eqn 34.2 & 34.3。
// 用于估计 SN shock cooling 的简单幂律函数。高于 10^(7.3) K 的部分是轫致辐射主导的。
Real Draine_2011::CoolingCurve(Real T_cgs) const {
  if (/*  1e5 < T_cgs &&  */ T_cgs <= std::pow(10.0, 7.3) ) {
    return 1.1e-22 * std::pow(T_cgs/1e6, -0.7) * std::exp(-1.18348e5/T_cgs); //TEMP 指数截断，具体的值有待调整
  } else if ( std::pow(10.0, 7.3) < T_cgs ) {  // 目前这里对 T 的上界没有限制，暂时没出问题
    return 2.3e-24 * std::pow(T_cgs/1e6, 0.5);
  }
  return 0.0;  // 默认返回值
}


Real Draine_2011_cutoff::CoolingCurve(Real T_cgs) const {
  if (T_cgs <= std::pow(10.0, 7.3) ) {
    return 1.1e-22 * std::pow(T_cgs/1e6, -0.7) * std::exp(- std::pow(4.5e4/T_cgs, 1.2)); // 这里指数截断的系数值是我自己手动拟合的。主要用于近似 Draine 2011 在 1e4 附近的截断
  } else if ( std::pow(10.0, 7.3) < T_cgs ) {  // 目前这里对 T 的上界没有限制，暂时没出问题
    return 2.3e-24 * std::pow(T_cgs/1e6, 0.5);
  }
  return 0.0;  // 默认返回值
}


Real OtherModel::CoolingCurve(Real T_cgs) const {
  throw std::invalid_argument("OtherModel is not implemented yet.");
}

Real OtherModel::Jacobian(Real T_cgs) const {
  throw std::invalid_argument("OtherModel is not implemented yet.");
}