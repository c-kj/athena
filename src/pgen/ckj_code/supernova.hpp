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
#include "../../parameter_input.hpp"

#include "region.hpp"


using Vector = std::array<Real, 3>;

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










#endif // SUPERNOVA_HPP_