#include <iostream>

#include "../../globals.hpp"
#include "../../task_list/task_list.hpp"

#include "ckj_plugin.hpp"

namespace ckj_plugin {
  // 这些全局变量需要定义
  const TimeIntegratorTaskList *ptlist = nullptr;  // const 防止修改指针指向的对象
  int current_stage = 0;
  int nstages = 0;
  Real source_term_weight;

  std::vector<Real> source_term_weights;  // 只在源文件中定义，不暴露给 namespace，免得调用时忘了 stage-1，导致越界 UB。

  void GetPointers(Outputs *pouts, TimeIntegratorTaskList *ptlist) {
    ckj_plugin::ptlist = ptlist;
    nstages = ptlist->nstages;

  }

  void UserWorkBeforeStage(Mesh *pmesh, int stage) {
    current_stage = stage;
    // 在主循环的每个 stage 中，都更新一下 source_term_weight，从而其他模块可以直接引用这个值。
    source_term_weight = source_term_weights[current_stage-1];  //* 注意 stage 从 1 开始，而 vector 从 0 开始，所以要 -1 。
  }

  // 计算指定 stage 的 source term weight。
  // 方法是：仿照 Athena++ 的 Integrator 对 cons 寄存器的计算步骤，只在所计算的 stage 注入大小为 1 的源项，从而最终获得的 cons 就是该 stage 的系数。
  // 计算步骤参见 TimeIntegratorTaskList::IntegrateHydro，省去对 Flux 和 坐标源项的处理。
  // 需要拿到 ptlist->stage_wghts，通过向源代码中的 TimeIntegratorTaskList 添加了一个函数 GetStageWeight 来获得。
  // 工具函数（只在源文件中定义，不暴露给 namespace）
  // 应该适用于各种 integrator，测试了 vl2, rk2, rk3, rk4 的无 orbital advection 情况。由于我也检查了 main_stage，也许对于 orbital advection 也适用，但不确定。
  // 同样不确定是否适用于具有特殊处理的 ssprk5_4，不过似乎 ssprk5_4 只是额外添加了 u2 寄存器，涉及 gamma_3 系数，但没有改变源项对 u 寄存器的作用。
  Real GetSourceTermWeight(int stage) {
    auto stage_wghts = ptlist->GetStageWeight();
    Real u = 0, u1 = 0;  // 初值都设为 0 ，从而看出源项的贡献
    for (int s = 1; s <= nstages; ++s) { // 模仿 Integrator 的计算步骤，走完所有的 stage，在每个 stage 中进行寄存器 u, u1 的运算。
      auto stage_wght = stage_wghts[s-1];
      if (stage_wght.main_stage) {  // 必须和 IntegrateHydro 一致，只在 main_stage 中进行寄存器运算！因为 nstages 中有一些不是 main_stage 的。
        u1 = u1 + stage_wght.delta * u;
        u  = stage_wght.gamma_1 * u + stage_wght.gamma_2 * u1 + (s==stage ? 1 : 0);   // 只在所求的 stage 注入大小为 1 的源项，从而结果就是该项的系数
      }
    }
    return u;
  }

  void GetSourceTermWeightList() {
    for (int stage = 1; stage <= nstages; ++stage) {
      source_term_weights.reserve(nstages);
      source_term_weights.push_back(GetSourceTermWeight(stage));
    }
    
    // 输出 source_term_weights 的信息
    if (Globals::my_rank == 0) {
      std::cout << "# ckj_plugin: calculating source_term_weights..." << '\n';
      std::cout << "nstages = " << nstages << '\n';
      std::cout << "Source term weights = ";
      for (auto w : source_term_weights) {
        std::cout << w << ", ";
      }
      std::cout << '\n';
    }
  }


}


// 这里用 __attribute__((weak)) 弱符号 提供默认实现。这样，如果 pgen 中没有实现，就会调用这里的默认实现，避免插入了 ckj_plugin 的 main.cpp 编译不通过。
// 目前不启用
#if false
namespace ckj_plugin {
    void __attribute__((weak)) GetPointers(Outputs *pouts, TimeIntegratorTaskList *ptlist) {
        // Do nothing
    }

    void __attribute__((weak)) UserWorkBeforeStage(Mesh *pmesh, int stage) {
        // Do nothing
    }

}
#endif