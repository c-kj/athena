#pragma once

#include <vector>
#include <algorithm>

#include "../../athena.hpp"

#include "cooling.hpp"

struct Cooling; // forward declaration

// 实现 Townsend 2009 算法的类
//* 这个类中所有的 T 和 Lambda 都是 cgs 单位制下的！就不单独标明了
struct TownsendCooling {
  TownsendCooling(Cooling *cooling, std::vector<Real> log10_T_array, std::vector<Real> log10_Lambda_array);
  // 委托构造函数，根据 CoolingFunc 在 log10_T_array 的采样点来构造
  // CoolingFunc 是一个接收 T_cgs，返回 Lambda_cgs （emissivity in cgs unit）的函数
  TownsendCooling(Cooling *cooling, std::vector<Real> log10_T_array, std::function<Real(Real)> CoolingFunc);

  Real TEF(Real T) const;  // TEF: Temporal Evolution Function
  Real TEF_inv(Real Y) const;

  int find_T_bin(Real log10_T) const;
  int find_Y_bin(Real Y) const;
  Real TEF_no_offset(Real log10_T, int i) const;

  // 计算新的温度 T_new。给定旧的温度 T (单位为 K)、密度 rho & 时间步长 dt (code unit)
  Real New_T(Real T, Real rho, Real dt) const;

  Cooling *pcooling; // Cooling 对象的指针，用于获取 mu_e、X_H 等参数
  Units *punit;  // Units 对象的指针


  const int N; // cooling curve 的点数。 = bin 数目 + 1

  std::vector<Real> log10_T_array, log10_Lambda_array; // Cooling curve 的数据，单位制为 cgs，以 log10 记录，为了计算更快，且与文献表格格式相同。长度为 N。
  std::vector<Real> alpha_array;  // 长度为 N。前 N-1 个元素是每个区间内的 slope，最后一个元素由 alpha_last 给定。
  std::vector<Real> Y_array;  // TEF 函数在各个节点的值。长度为 N。
  std::vector<Real> T_Lambda_factor_array;  // 用于 TEF 中的因子 (T_k / T_N) / (Lambda_k / Lambda_N)，预先算好。长度为 N。

  const Real alpha_last = 0.5;  // 最后一个区间的 alpha，默认为 0.5 （free-free 轫致辐射）
  Real log10_T_ref, log10_Lambda_ref;  // 参考点
  Real T_min; // cooling curve 数据中的最低温度，cgs 单位制

  Real rho_dt_coeff;  // New_T 函数中乘在 rho * dt 前面的系数

private:
  // 辅助函数，根据 CoolingFunc 计算对应于每个 log10_T 的 log10_Lambda，同时检查是否有 Lambda <= 0 的情况
  // CoolingFunc 是一个接收 T_cgs，返回 Lambda_cgs 的函数
  std::vector<Real> Calc_log10_Lambda_array(std::function<Real(Real)> CoolingFunc, std::vector<Real> log10_T_array);
};


// 寻找 x 属于哪个 bin，返回 bin 的 index（左边界在 bin_edges 中的 index）。
// bin_edges 是 bin 的边界，每个 bin 视为前闭后开区间 [ bin_edges[i], bin_edges[i+1] )。注意若 cmp 为 > ，则数值大的那一侧是闭的。
// 假定：bin_edges 在 cmp 下是有序的，且 size >= 2，否则 UB。默认的 cmp 是 std::less<T>()，即从小到大排序。
// 边界情况处理：
// 如果 x 的值在最后一个 bin 的右边界之外，不做特殊处理，正常返回最后一个 bin_edge 的 index (视为处在最后一个延伸到无穷的 bin 中)。
// 如果 x 的值在第一个 bin 的左边界之外，返回 -1，由外部处理。
template<typename T, typename Compare = std::less<T>>
int find_bin_index(const std::vector<T>& bin_edges, const T& x, Compare cmp = Compare()) {  // cmp 的默认值是 std::less<T>()
  // 使用「前后」或「左右」来描述，而非「大小」，因为 cmp 可能是任意比较函数
  auto it_after = std::upper_bound(bin_edges.begin(), bin_edges.end(), x, cmp); // 指向 bin 区间的右边界
  auto it_before = std::prev(it_after); // 指向 bin 区间的左边界

  int i = std::distance(bin_edges.begin(), it_before);
  return i;
}