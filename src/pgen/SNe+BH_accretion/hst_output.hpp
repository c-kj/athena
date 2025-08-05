#pragma once

#include <memory> // std::unique_ptr

#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"

/* -------------------------------------------------------------------------- */
/*                               自定义的 hst 机制                              */
/* -------------------------------------------------------------------------- */

//* 使用方法：
// 1. 在 Mesh::InitUserMeshData 中调用 HSTManager 的构造函数，为全局变量 hst 赋值。（构造时会自动 Enroll 所有的 hst 条目）
// 2. 在 MeshBlock::InitUserMeshBlockData 中 AllocateRealUserMeshBlockDataField 后，调用 hst->RequestMeshBlockData(this) 来初始化 hst 机制所使用的 MB_data.
// 3. 在 pgen 中定义 HSTManager::register_all_entries 和 HSTManager::UserWorkBeforeHstOutput 接口
// 然后添加自定义的 hst 条目，步骤如下：

//* 添加新的 hst 条目的步骤
// 1. 在 HSTManager::register_all_entries 中，使用 add_entry 添加自定义的 hst output 条目
// 2. 添加相应的 hst 输出所需要的计算步骤 (对 hst->get_proxy(pmb) 进行操作)
//  - 瞬时值：在 HSTManager::UserWorkBeforeHstOutput 直接赋值
//  - 瞬时统计量：在 HSTManager::UserWorkBeforeHstOutput 中计算，赋值
//  - 累积量：在相应的源项中进行累计。别忘了乘以相应 stage 的 source_weight。


// Forward declaration
struct HSTManager;


// 表示一个 hst output 的元信息条目
struct HSTEntry {
  std::string name;
  UserHistoryOperation op;
};

// 代理 HSTManager 的 MB_data，[] 运算符接收 hst 的 name 而非 index，从而更方便访问、修改
struct HSTProxy {
  HSTManager *phst;  // 所属的 HSTManager
  AthenaArray<Real>& data;  // 对 HSTManager 的 MB_data 的引用
  
  HSTProxy(HSTManager *phst, AthenaArray<Real>& data);

  Real& operator[](const std::string& name);
};

/*
提供一个通用的 hst 机制框架。特定于具体 pgen 的内容，由 pgen 中定义相应的接口。
- 简化了原有的 Enroll 过程。只需输入 name 和 op（默认为 sum），不需要手动输入 index，避免错乱；也不需要单独定义函数，避免繁冗。
- 通过 HSTProxy 实现以字符串 name 作为索引，而非整数 index，更直观不易错
  - 可以使用程序生成字符串，从而程序化地组合、批量定义/赋值，大大提高可扩展性
- 利用 ruser_meshblock_data 来储存 hst output，对于累积量，还能自动跨 MPI 汇总、跨 restart 恢复。
  - 所有的数据储存在以 name 为索引的一维数组中，不再需要单独建立索引再复制过来。
  - 自动注册、引用相应的 ruser_meshblock_data
- 使用 UserWorkBeforeHstOutput，只在需要计算 hst output 时调用。并且在一次遍历中计算所有瞬时统计量，避免重复进入循环。
*/
struct HSTManager {
  
  /* ---------------------------------- 用于初始化 --------------------------------- */
  
  // 构造函数。在 Mesh::InitUserMeshData 中调用。
  HSTManager(Mesh *pmesh, int MB_data_index);
  
  // 在 MeshBlock::InitUserMeshBlockData 中调用
  // 用于在每个 MeshBlock 上开辟 HSTManager 机制所需的 ruser_meshblock_data 数组的内存
  void RequestMeshBlockData(MeshBlock *pmb) {
    MB_data(pmb).NewAthenaArray(size());
  }
  
  /* ---------------------------------- 公开接口 ---------------------------------- */
  
  // 返回 hst output 条目的数量
  int size() const { return static_cast<int>(entries.size()); }
  
  // HSTManager 机制在每个 MeshBlock 上所使用的数组 ruser_meshblock_data 的引用
  // 此接口用于便捷、统一地访问、修改这个数组
  // 不过其实在 pgen 中不太需要直接访问这个数组，而是通过 HSTProxy 来访问。也可改为 private。
  AthenaArray<Real>& MB_data(MeshBlock *pmb) {
    return pmb->ruser_meshblock_data[MB_data_index];
  }
  
  // 创建一个 proxy，从而便捷地用 name 访问 MB_data
  HSTProxy get_proxy(MeshBlock *pmb) {
    return HSTProxy(this, MB_data(pmb));
  }
  
  /* ---------------------------------- 成员变量 ---------------------------------- */
  
  const int MB_data_index;          // hst 机制所使用的数组 在 ruser_meshblock_data 中的 index
  std::vector<HSTEntry> entries;    // 添加到此 HSTManager 的所有 hst output 条目，按添加顺序排列
  std::unordered_map<std::string, int> name_index_map;  // 用于根据 name 查找对应的 index
  
private:  // 防止外部调用
  
  /* ----------------------------- 必须由 pgen 定义的接口 ----------------------------- */
  // 只能在 pgen 中定义，但调用时，要用公开的接口（包装了以下接口）
  
  // 在这个函数中，使用 add_entry 添加所有自定义的 hst output 条目
  void register_all_entries(); 
  
  // 在每个 hst 输出前调用。
  // 可用于计算瞬时量、瞬时统计量（避免重复进入循环）等
  void UserWorkBeforeHstOutput(MeshBlock *pmb);
  
  /* --------------------------------- 私有成员函数 --------------------------------- */
  
  // 添加新的 hst output 条目。只应该在 HSTManager::register_all_entries 中调用
  void add_entry(std::string name, UserHistoryOperation op=UserHistoryOperation::sum) {
    entries.emplace_back(HSTEntry{std::move(name), op});
  }
  
  // 在构造函数中调用，把所有的 hst output enroll 到 pmesh 中。
  // 注意：重复调用这个函数会重新 allocate，清空之前 enroll 的 hst 信息
  void enroll_all_entries(Mesh *pmesh);
  
  // 统一用这个函数，根据不同的 iout，给出不同的 hst 的数据。
  // 这个函数必须是 static 的，不能是成员函数，否则无法传入给 Mesh::EnrollUserHistoryOutput
  // 另外，还在每个 cycle 中调用 UserWorkBeforeHstOutput
  static Real hst_function(MeshBlock *pmb, int iout);  
};

// 全局的 hst 机制指针。在 hst_output.cpp 中定义。其他文件可以直接使用这个指针。
// 需要在 Mesh::InitUserMeshData 中赋值初始化。
// 这个类的设计上应该是单例（在每个 MPI rank 上），但并不做任何检查来保证。
//* 注意：由于 hst_function 中调用了 hst，所以其他地方都必须使用「hst」这个指针。
extern std::unique_ptr<HSTManager> hst;