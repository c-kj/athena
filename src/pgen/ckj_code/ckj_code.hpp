#pragma once

// 这里放一些用于在多个自定义文件中 include 的声明。

#include <map>

// Passive Scalar 与其对应的 index 的对应
enum PassiveScalarIndex {SN};

// SourceTerm （主要是 SN 和 cooling）的注入时机。按照单个 cycle 内的先后顺序排列，从而可以进行 < 比较
// 有关源项的时机，参看我的笔记： 《黑洞吸积 SN 湍流 project.md》和《黑洞吸积 SN 湍流 TODO 历史.md》
enum class SourceTermPosition {InSourceTerm, AfterSourceTerm, UserWorkInLoop}; 

// 用于处理 input 参数中的字符串。对这个 map 使用 .at 方法来获取对应的 SourceTermPosition，如果不存在则会抛出异常
//? 这里也许不能声明为 static?
static std::unordered_map<std::string, SourceTermPosition> source_term_position_map = {
  {"InSourceTerm", SourceTermPosition::InSourceTerm},
  {"AfterSourceTerm", SourceTermPosition::AfterSourceTerm},
  {"UserWorkInLoop", SourceTermPosition::UserWorkInLoop}
};

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