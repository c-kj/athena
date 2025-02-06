
#include "Townsend.hpp"
#include "cooling.hpp"

TownsendCooling::TownsendCooling(Cooling *cooling, std::vector<Real> log10_T_array, std::vector<Real> log10_Lambda_array): 
    pcooling(cooling), punit(cooling->punit), 
    N(log10_T_array.size()), log10_T_array(log10_T_array), log10_Lambda_array(log10_Lambda_array) {
  // 检查形状
  if (N != log10_Lambda_array.size()) {
    throw std::runtime_error("log10_T_array 和 log10_Lambda_array 的长度不一致！");
  }
  // 检查 log10_T_array 是否递增
  if (not std::is_sorted(log10_T_array.begin(), log10_T_array.end())) {
    throw std::runtime_error("log10_T_array 不是递增的！");
  }
  
  T_min = std::pow(10, log10_T_array[0]);

  //* Townsend 2009 文献中的下标 k 是从 1 到 N-1 的，而 C++ 中的 i 是从 0 到 N-2 的

  // 初始化 alpha_array
  alpha_array.resize(N);
  for (int i = 0; i < N-1; ++i) {
    alpha_array[i] = (log10_Lambda_array[i+1] - log10_Lambda_array[i]) / (log10_T_array[i+1] - log10_T_array[i]);
  }
  alpha_array[N-1] = alpha_last;

  // 取最大的 T 作为参考点
  log10_T_ref = log10_T_array[N-1];  
  log10_Lambda_ref = log10_Lambda_array[N-1];  

  // 初始化 T_Lambda_factor_array
  T_Lambda_factor_array.resize(N);
  for (int i = 0; i < N; ++i) {
    T_Lambda_factor_array[i] = std::pow(10, (log10_T_array[i] - log10_T_ref) - (log10_Lambda_array[i] - log10_Lambda_ref) );
  }

  // 初始化 Y_array
  Y_array.resize(N);
  Y_array[N-1] = 0.0;
  for (int i = N-2; i >= 0; --i) {
    Y_array[i] = Y_array[i+1] - TEF_no_offset(log10_T_array[i+1], i);
  }
  // 检查 Y_array 是否单调递减
  if (not std::is_sorted(Y_array.begin(), Y_array.end(), std::greater<Real>())) {
    throw std::runtime_error("Y_array 不是递减的！");
  }


  // 计算 New_T 中所需的系数 rho_dt_coeff
  Real gamma = pcooling->gamma;
  Real Lambda_ref = std::pow(10, log10_Lambda_ref);
  Real T_ref = std::pow(10, log10_T_ref);
  Real k_B_cgs = Constants::k_boltzmann_cgs;
  Real m_H_cgs = Constants::hydrogen_mass_cgs;
  constexpr Real mu = Abundance::mu;  // 目前这些元素丰度都是常量。以后可能再考虑空间非均匀的情况
  constexpr Real X_H = Abundance::X_H;
  constexpr Real mu_e = Abundance::mu_e;
  rho_dt_coeff = Lambda_ref/T_ref * (gamma-1) / (k_B_cgs * m_H_cgs) * (mu * X_H / mu_e); // 系数的表达式，自己手推的。这里所有量都是 cgs 单位下的。
  rho_dt_coeff *= punit->code_density_cgs * punit->code_time_cgs; // 提前乘上 rho*dt 的单位转换系数
}


TownsendCooling::TownsendCooling(Cooling *cooling, std::vector<Real> log10_T_array, std::function<Real(Real)> CoolingFunc):
  TownsendCooling(cooling, log10_T_array, Calc_log10_Lambda_array(CoolingFunc, log10_T_array)) {}


Real TownsendCooling::TEF(Real T) const {
  Real log10_T = std::log10(T);
  int i = find_T_bin(log10_T);
  // 这里不对 i 进行检查。在外部已经处理了 T < T_min 的情况；而如果 T > T_max，那么会落入 [T_max, \infty) 的区间，使用 alpha_array[N-1] = 0.5 的轫致辐射。
  return Y_array[i] + TEF_no_offset(log10_T, i);
}

Real TownsendCooling::TEF_inv(Real Y) const {
  int i = find_Y_bin(Y);
  if (i == -1) {return T_min;}  // 如果 Y > Y_array[0] 越界，说明对于给定的 Y（约化 dt），在 T range 内找不到对应的 T，达到了 cooling curve 的最低温度。返回 T_min

  Real alpha_k = alpha_array[i];
  Real log10_T_k = log10_T_array[i];

  Real Y_diff = Y - Y_array[i];

  Real T_Lambda_factor = T_Lambda_factor_array[i];
  Real T_Lambda_factor_inv = 1 / T_Lambda_factor;

  return std::pow(10, log10_T_k) * (alpha_k == 1 ? 
    std::exp(- T_Lambda_factor_inv * Y_diff) : 
    std::pow(1 - (1 - alpha_k) * T_Lambda_factor_inv * Y_diff, 1 / (1 - alpha_k)) 
  );
}



// 幂函数的倒数的定积分 \int_T^{T_ref} T^{-\alpha} dT
// 这里 ln_T_diff 是 ln(T / T_k)。初始化 Y_k 时， T 取 T_{k+1}。一般情况下 ln_T_diff > 0。
inline Real powerlaw_integral(Real ln_T_diff, Real alpha_k) {
  return alpha_k == 1 ? -ln_T_diff : 1 / (1 - alpha_k) * (1 - std::exp(ln_T_diff * (1 - alpha_k)) );
}

// 输入 T（取 log10）和指定的 bin，计算不含 Y_k 的 TEF（Townsend 2009 公式 A5 的大括号部分）
// 这里的 i 起到公式中 k 的作用（但起点不同），标记当前 bin。i 不能自动从 T 推断，因为在初始化 Y_k 时，输入的 T 是当前 bin 的上端点。
// 用于复用：在 TEF 中使用时，加上 Y_k，并且 i 通过 find_T_bin 计算得到；在初始化各个 bin 中的 Y_k 时，。
Real TownsendCooling::TEF_no_offset(Real log10_T, int i) const {
  // 这里不对 i 进行检查，避免影响性能。在外部已经检查了。
  Real alpha_k = alpha_array[i];
  Real log10_T_k = log10_T_array[i];

  Real ln_T_diff = (log10_T - log10_T_k) * std::log(10);
  Real T_Lambda_factor = T_Lambda_factor_array[i];  

  return T_Lambda_factor * powerlaw_integral(ln_T_diff, alpha_k);
}


//TODO 对于等间距的情况，直接计算 index 更快，不用二分查找
int TownsendCooling::find_T_bin(Real log10_T) const {
  return find_bin_index(log10_T_array, log10_T);
}

int TownsendCooling::find_Y_bin(Real Y) const {
  return find_bin_index(Y_array, Y, std::greater<Real>()); // Y_array 是递减的
}

Real TownsendCooling::New_T(Real T, Real rho, Real dt) const {
  if (T <= T_min) {return T_min;}  // 如果输入的 T 已经低于 cooling curve 数据的最低温度，直接返回最低温度
  return TEF_inv(TEF(T) + rho_dt_coeff * rho * dt); // 这里的 rho*dt 是 code unit 下的，但 rho_dt_coeff 中已经乘了对应单位的转换系数，其整体是在 cgs 下进行的计算。
}

std::vector<Real> TownsendCooling::Calc_log10_Lambda_array(std::function<Real(Real)> CoolingFunc, std::vector<Real> log10_T_array) {
  std::vector<Real> log10_Lambda_array;
  for (Real log10_T: log10_T_array) {
    Real T = std::pow(10, log10_T);
    Real Lambda = CoolingFunc(T);
    if (Lambda <= 0) {throw std::invalid_argument("Cooling rate Lambda 必须 >= 0 !");}  // 不支持 Lambda == 0 的情况
    log10_Lambda_array.push_back(std::log10(Lambda));
  }
  return log10_Lambda_array;
}