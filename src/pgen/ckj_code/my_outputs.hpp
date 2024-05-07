#ifndef MY_OUTPUTS_HPP_
#define MY_OUTPUTS_HPP_

// uov 的索引。所有使用 uov 数组的地方都必须使用这里的 enum 值而非自己填 int
enum UserOutputVariableIndex {
  // I_test,
  I_cooling_rate,
  // 其他 uov ...
  N_UOV  // uov 的个数，用于 AllocateUserOutputVariables。必须放在 enum 的最后
};

#endif  // MY_OUTPUTS_HPP_