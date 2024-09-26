#ifndef CKJ_PLUGIN_HPP
#define CKJ_PLUGIN_HPP

#pragma once

// 目前没有做 forward declaration，但似乎还能正常编译，因为每个 include 这个头文件的 .cpp 都会包含所需？
class Outputs;
class TaskList;
class Mesh;


// 自定义的 plugin 函数，（目前）用于插入到 main.cpp 中，获取一些拿不到的变量，比如指针、stage 等。
// 这些函数应当在自己的 pgen 文件中实现。
// 如果 pgen 中没有实现，则需要在 ckj_plugin.cpp 中，用弱符号 __attribute__((weak)) 实现一个默认版本。
namespace ckj_plugin { 

  // 获取一些拿不到的指针，用于我自定义的模块中的操作。
  // 调用时机：在 main.cpp 中，主循环开始前。
  // 注意：在主循环开始前，比如 Mesh 初始化时、ProblemGenerator 设定初始条件时，这个函数还没有获取到指针！
  void GetPointers(Outputs *pouts, TaskList *ptlist);

  // 调用时机：在每个 cycle 的每个 stage 开始时
  void UserWorkBeforeStage(Mesh *pmesh, int stage);
  
  // 后续可以扩充，比如拿到 Outputs 的指针、TaskList 的指针等等。
}

#endif // CKJ_PLUGIN_HPP