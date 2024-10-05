#ifndef CKJ_PLUGIN_HPP
#define CKJ_PLUGIN_HPP

#pragma once

#include <vector>

#include "../../athena.hpp"

// 目前没有做 forward declaration，但似乎还能正常编译，因为每个 include 这个头文件的 .cpp 都会包含所需？
class Outputs;
class TimeIntegratorTaskList;
class Mesh;


// 自定义的 plugin 函数，（目前）用于插入到 main.cpp 中，获取一些拿不到的变量，比如指针、stage 等。
// 这些函数应当在自己的 pgen 文件中实现。
// 如果 pgen 中没有实现，则需要在 ckj_plugin.cpp 中，用弱符号 __attribute__((weak)) 实现一个默认版本。
// 单独开一个命名空间，是为了在引用时更明确地看出这是 plugin，而非不知哪里冒出来的全局变量。
namespace ckj_plugin { 

  extern const TimeIntegratorTaskList *ptlist;   // 从 main.cpp 中获取的指针
  extern int current_stage;                // 当前的 stage，从 main.cpp 的 for 循环中获取
  extern int nstages;                      // 所选的 time integrator 的 stage 总数，用于判断是否为最后一个 stage。从 ckj_plugin::GetPointers 中获取。
  extern Real source_term_weight;          // 当前 stage 对应的源项权重，即该 stage 内源项函数对 cons 造成的改变量最终进入 cycle 结尾的 cons 的比例。注意：这些「weight」之和往往不等于 1。

  // 获取一些拿不到的指针，用于我自定义的模块中的操作。
  // 调用时机：在 main.cpp 中，主循环开始前。
  // 注意：在主循环开始前，比如 Mesh 初始化时、ProblemGenerator 设定初始条件时，这个函数还没有获取到指针！
  void GetPointers(Outputs *pouts, TimeIntegratorTaskList *ptlist);

  // 调用时机：在每个 cycle 的每个 stage 开始时
  void UserWorkBeforeStage(Mesh *pmesh, int stage);

  // 用于初始化 source_term_weights
  // 调用时机：在 main.cpp 中，主循环开始前。这假定了 integrator 在主循环中不变。不确定 restart 时更改会怎样。
  void GetSourceTermWeightList();
}

#endif // CKJ_PLUGIN_HPP