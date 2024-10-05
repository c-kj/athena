#ifndef MY_OUTPUTS_HPP_
#define MY_OUTPUTS_HPP_

//* 注意 enum 的枚举成员是全局的，要避免命名冲突！如果不好避免，则用 namespace 包裹

// uov 的索引。所有使用 uov 数组的地方都必须使用这里的 enum 值而非自己填 int
namespace UOV {
  enum UserOutputVariableIndex {
    // test,
    cooling_rate,
    // 其他 uov ...
    N_UOV  // uov 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
  };
}

// 名字很长，但这样比较明确。引用时可以临时给个别名，比如 namespace idx = RealUserMeshBlockDataIndex;
namespace RealUserMeshBlockDataIndex {
  enum RealUserMeshBlockDataIndex {
    // 源项造成的改变量的历史记录
    // 吸积量
    accreted_mass,
    accretion_rate,
    accreted_SN_tracer,
    accretion_rate_SN_tracer,
    accreted_energy,
    accreted_momentum,
    accreted_angular_momentum,
    // 能量
    total_cooling_loss,
    BH_gravity_work,
    // SN 注入量
    SN_injected_energy,
    SN_injected_mass,
    SN_injected_number,
    // 其他 ...
    N_RealUserMeshBlockData  // DataField 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
  };
}

#endif  // MY_OUTPUTS_HPP_