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



// forward declaration
struct Supernovae;
struct SupernovaParameters;
struct SupernovaEvent;


// 表示单个 SN 爆炸的事件。只包含简单的信息（t,x,v 等），其他统一的参数在 SupernovaParameters 对象中。 
struct SupernovaEvent {
  // 构造函数，这里只是简单地初始化成员变量
  SupernovaEvent(SupernovaParameters *paras, Real t, Point x, Vector v); 

  SupernovaParameters *paras;  // 该事件所属的 SupernovaParameters
  Real time;
  Point position;
  Vector velocity;
  Real velocity_magnitude;

  // energy 和 mass 注入的 Region。
  // 之所以用指针，是因为这样可以适用于各种 Region 类型，不需要指定具体是哪个子类。用智能指针，不需要手动释放内存
  std::unique_ptr<Region> energy_region = nullptr; 
  std::unique_ptr<Region> mass_region = nullptr;
  Real energy_density = 0.0, mass_density = 0.0;

  std::string Info();
};


// 单个（种）超新星的信息。目前，一个 SupernovaParameters 对象代表「一种超新星」，即所有参数相同、只有时间和位置不同的多个超新星。
struct SupernovaParameters {
  SupernovaParameters(Supernovae *pSNe, ParameterInput *pin, const int i);
  
  const int ndim;
  // 基本类型 (int, double, ...) 的成员变量需要显式初始化，否则是未定义的值。
  Real energy = 0.0;
  Real mass = 0.0;

  std::string region_name;  //TODO 名称改为 injection_shape ？连着输入参数模板一起改。 反正输入参数中目前还无需指定，只用默认值。
  Real radius, mass_radius;

  std::vector<std::unique_ptr<SupernovaEvent>> event_list;

  // 不一定会用上的属性。目前还没写实现
  // std::string name;
  // std::string type;
    
};



// 用于存储所有超新星，控制其注入
struct Supernovae {
  Supernovae() = default; // 默认构造函数。需要这个才能在声明 Supernovae supernovae 时初始化。不过如果改用指针，可能就不需要这个了。
  Supernovae(Mesh *pmy_mesh, ParameterInput *pin);

  Mesh *pmy_mesh;
  Units *punit;  // TEMP 暂时放在这个类中，考虑挪到 SupernovaParameters 中
  int ndim;      // TEMP 暂时放在这个类中，考虑挪到 SupernovaParameters 中
  int SN_flag;
  
  // 目前直接把这个量设为编译时常量，不从参数文件读取。SN 在 SourceTerm 结尾注入时最佳的，能够准确控制时间并且能在当前步的 prim 中输出。
  static constexpr SourceTermPosition source_term_position = SourceTermPosition::AfterSourceTerm; // 在哪个位置注入 SN

  //TODO 挪到 SupernovaParameters 类中
  Real SN_energy_unit;
  Real SN_mass_unit;

  std::vector<std::unique_ptr<SupernovaParameters>> supernova_paras_list;  // 储存所有 SupernovaParameters。这里用 unique_ptr，表明所有权在此管理。
  // 以下两个 vector 使用 raw pointer，因为 SupernovaEvent 由 SupernovaParameters 管理。
  std::vector<SupernovaEvent*> supernova_list;  // 储存所有单个 SupernovaEvent
  std::vector<SupernovaEvent*> supernova_to_inject;  // 储存当前时刻需要注入的超新星，在每个主循环的 cycle 中更新

  void InitSupernovaParameters(ParameterInput *pin);
  void InitSupernovaEvents();
  void PrintInfo();

  void GetSupernovaeToInject(const Real time, const Real dt);
  void SuperNovaeSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

  Real SupernovaeTimeStep(MeshBlock *pmb);

};







#endif // SUPERNOVA_HPP_