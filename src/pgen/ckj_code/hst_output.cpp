#include "hst_output.hpp"

HSTProxy::HSTProxy(HSTManager *phst, AthenaArray<Real>& data): 
  phst(phst), data(data) {}

Real& HSTProxy::operator[](const std::string& name) {
  // 对传入的 name 进行检查，如果不存在，直接报错退出。用于检测 拼写错误、忘记 add_entry、条件不对 等 bug。
  // 这里不使用异常而是直接 abort，因为这个错误不应当被捕获。这样对性能影响较小。
  if (phst->name_index_map.find(name) == phst->name_index_map.end()) {
      std::cerr << "HSTProxy: name not found: " << name << '\n';
      std::abort();
  }

  return data(phst->name_index_map[name]);
}



HSTManager::HSTManager(Mesh *pmesh, int MB_data_index): MB_data_index(MB_data_index) {
  register_all_entries();  // 添加所有自定义的 hst output 条目

  // 建立 name_index_map
  for (int i = 0; i < entries.size(); ++i) {
      const auto& entry = entries[i];
      name_index_map[entry.name] = i;
  }

  enroll_all_entries(pmesh);  // 把所有 hst output 条目注册到 pmesh 上
}


// 由于需要调用 Mesh 的 private 成员函数，所以需要在 Mesh 类的声明中加一行 friend class HSTManager;
// 如果不修改 Mesh 类，则需要把这个函数的内容展开到 Mesh::InitUserMeshData 中（注意把 this 改为 hst, pmesh 改为 this）
void HSTManager::enroll_all_entries(Mesh *pmesh) {
  pmesh->AllocateUserHistoryOutput(size());
  for (int i = 0; i < size(); ++i) {
    const auto& entry = entries[i];
    pmesh->EnrollUserHistoryOutput(i, hst_function, entry.name.c_str(), entry.op);
  }
}


Real HSTManager::hst_function(MeshBlock *pmb, int iout) {
  //* 这里的 hst 是那个单例的全局变量。因为这个函数不能是一个成员函数，所以只能用这种方式访问。
  // 每个 MeshBlock 的每个 cycle 中，第一次要输出 hst output 时，调用 UserWorkBeforeHstOutput
  //* 注意：如果没有添加任何一个自定义的 hst output，就不会调用这个函数，从而不会执行 hst->UserWorkBeforeHstOutput
  if (iout == 0) {hst->UserWorkBeforeHstOutput(pmb);} 
  // 根据 iout 返回不同的 hst 数据
  return hst->MB_data(pmb)(iout);
}

std::unique_ptr<HSTManager> hst;