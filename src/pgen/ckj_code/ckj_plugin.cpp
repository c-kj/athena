#include "ckj_plugin.hpp"
#include "../../task_list/task_list.hpp"


namespace ckj_plugin {
  // 这些全局变量需要定义
  TimeIntegratorTaskList *ptlist = nullptr;
  int current_stage = 0;
  int nstages = 0;

  void GetPointers(Outputs *pouts, TimeIntegratorTaskList *ptlist) {
    ckj_plugin::ptlist = ptlist; // 规避重名，用 :: 来指代全局命名空间。
    nstages = ptlist->nstages;

  }

  void UserWorkBeforeStage(Mesh *pmesh, int stage) {
    current_stage = stage;
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