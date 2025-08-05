#ifndef MY_OUTPUTS_HPP_
#define MY_OUTPUTS_HPP_

//* 注意 enum 的枚举成员是全局的，要避免命名冲突！如果不好避免，则用 namespace 包裹

// uov 的索引。所有使用 uov 数组的地方都必须使用这里的 enum 值而非自己填 int
namespace UOV {
  enum UserOutputVariableIndex {
    // test,
    cooling_rate,
    heating_rate,
    // 其他 uov ...
    N_UOV  // uov 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
  };
}

// 名字很长，但这样比较明确。引用时可以临时给个别名，比如 namespace idx = RealUserMeshBlockDataIndex;
namespace RealUserMeshBlockDataIndex {
  enum RealUserMeshBlockDataIndex {
    hst,  // 用于 hst 机制
    // 其他 ...
    N_RealUserMeshBlockData  // DataField 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
  };
}


#endif  // MY_OUTPUTS_HPP_