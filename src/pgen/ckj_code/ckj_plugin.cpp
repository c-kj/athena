#include "ckj_plugin.hpp"


// 这里用 __attribute__((weak)) 弱符号 提供默认实现。这样，如果 pgen 中没有实现，就会调用这里的默认实现，避免插入了 ckj_plugin 的 main.cpp 编译不通过。
namespace ckj_plugin {
    void __attribute__((weak)) GetPointers(Outputs *pouts, TaskList *ptlist) {
        // Do nothing
    }

    void __attribute__((weak)) UserWorkBeforeStage(Mesh *pmesh, int stage) {
        // Do nothing
    }

}
