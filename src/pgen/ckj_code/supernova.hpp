// 实现超新星爆炸

#ifndef SUPERNOVA_HPP_
#define SUPERNOVA_HPP_


// C++ headers
#include <array>
#include <vector>
#include <string>
#include <memory>   // std::unique_ptr

// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"

// 自定义的头文件
#include "ckj_code.hpp"
#include "region.hpp"


using Vector = std::array<Real, 3>;

// 单个（组）超新星的信息
struct SuperNova {
  SuperNova(ParameterInput *pin, const int i, const int dim);
  
  std::vector<Real> time_list;
  Real energy = 0.0;
  Real mass = 0.0;
  Vector velocity = {0, 0, 0};

  // energy 和 mass 注入的 Region。
  // 之所以用指针，是因为这样可以适用于各种 Region 类型，不需要指定具体是哪个子类。用智能指针，不需要手动释放内存
  std::unique_ptr<Region> energy_region = nullptr; 
  std::unique_ptr<Region> mass_region = nullptr;

  std::string region_name;
  Real energy_density=0.0, mass_density=0.0;

  // 不一定会用上的属性。目前还没写实现
  // std::string name;
  // std::string type;
    
};


std::vector<SuperNova> read_supernova_list(ParameterInput *pin, const int ndim) ;


// 用于存储所有超新星，控制其注入
struct SuperNovae {
public:
  SuperNovae() = default; // 默认构造函数。需要这个才能在声明 SuperNovae supernovae 时初始化。不过如果改用指针，可能就不需要这个了。
  SuperNovae(ParameterInput *pin, Mesh *pmy_mesh);

  Units *punit;
  int ndim;
  int SN_flag;
  std::string integrator;
  
  // 目前直接把这个量设为编译时常量，不从参数文件读取。SN 在 SourceTerm 结尾注入时最佳的，能够准确控制时间并且能在当前步的 prim 中输出。
  static constexpr SourceTermPosition source_term_position = SourceTermPosition::AfterSourceTerm; // 在哪个位置注入 SN

  Real SN_energy_unit;
  Real SN_mass_unit;

  std::vector<SuperNova> supernova_list;
  std::vector<SuperNova*> supernova_to_inject;

  void GetSupernovaeToInject(const Real time, const Real dt);
  void SuperNovaeSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);


};







#endif // SUPERNOVA_HPP_