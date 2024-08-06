#ifndef MY_OUTPUTS_HPP_
#define MY_OUTPUTS_HPP_

//* 注意 enum 的枚举成员是全局的，要避免命名冲突！如果不好避免，则用 namespace 包裹

// uov 的索引。所有使用 uov 数组的地方都必须使用这里的 enum 值而非自己填 int
enum UserOutputVariableIndex {
  // I_test,
  I_cooling_rate,
  // 其他 uov ...
  N_UOV  // uov 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
};

enum RealUserMeshBlockDataIndex {
  I_accreted_mass,
  I_accretion_rate,
  I_accreted_SN_tracer,
  I_accretion_rate_SN_tracer,
  I_accreted_momentum,
  I_accreted_angular_momentum,
  // 其他 ...
  N_RealUserMeshBlockData  // DataField 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
};

#endif  // MY_OUTPUTS_HPP_