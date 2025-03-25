// 这个头文件用于在各个特定于 SNe+BH_accretion 这个 pgen 的模块（文件）之间，共享函数、全局变量的声明
// 如果只在 pgen 单个文件内使用，放在 cpp 文件中即可。如果要适用于多个 pgen，则放在 ckj_code.hpp 中。

#pragma once


/* -------------------------------------------------------------------------- */
/*                                   静态参数                                   */
/* -------------------------------------------------------------------------- */

//TODO 把 Passive Scalar Index 改为类似自定义 hst 机制那样，用字符串做索引，更灵活，可以适应 key 不存在的情况。另外可以写一个机制在编译时保证字符串索引不会出错。
// Passive Scalar 与其对应的 index 的对应
//* 修改时，记得在 pgen 中更改对应的名字
namespace PassiveScalarIndex {
  enum PassiveScalarIndex {
    SN,
    initial_radius,
    initial_fluid,
    N_PassiveScalar_defined
  };
}

static_assert(PassiveScalarIndex::N_PassiveScalar_defined == NSCALARS, "### ERROR: 已定义的 Passive Scalar 数目与 NSCALARS 不一致！");

// 元素丰度，用于计算 mu，从而计算温度 T。目前用在设定初值、计算 Cooling 上。
// 目前暂时都 hard code 为常量。
  // 如果要在 input file 中设定，则需要写构造函数，接收 ParameterInput *pin。
  // 如果要改为动态演化，那么就是逐点的了，需要重构。
struct Abundance {
  //* 这里修改的话，要在 python 后处理中同步修改。
  // 设定元素丰度（H, He, Metal 的质量分数）
  static constexpr Real X_H = 0.7, X_Metal = 0.01295;  //? 这是什么丰度？太阳的？
  // X_H = 1.0, X_Metal = 0.0;   // 纯 H
  static constexpr Real X_He = 1.0 - X_H - X_Metal;

  static constexpr Real mu = 4.0/(5*X_H + 3 - X_Metal);  //? mu 取什么值？
  static constexpr Real mu_e = 2.0/(1+X_H);
};


/* -------------------------------------------------------------------------- */
/*                                     函数                                     */
/* -------------------------------------------------------------------------- */
