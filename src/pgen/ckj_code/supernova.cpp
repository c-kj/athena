
#include <array>
#include <vector>
#include <sstream>  // std::stringstream, std::getline
#include <memory>   // std::unique_ptr
#include <stdexcept>  // std::invalid_argument

#include "supernova.hpp"
#include "utils.hpp"



SuperNova::SuperNova(ParameterInput *pin, const int i, const int dim) {
  std::string block_name, i_str, center_str;
  Real radius, energy_radius, mass_radius;
  Point center, energy_center, mass_center;
  block_name = "supernova";

  i_str = std::to_string(i);

  time_list = read_vector(pin->GetString(block_name, "SN_" + i_str + "_time"));

  energy = pin->GetOrAddReal(block_name, "SN_" + i_str + "_energy", 0.0);
  mass = pin->GetOrAddReal(block_name, "SN_" + i_str + "_mass", 0.0);
  velocity = read_array(pin->GetOrAddString(block_name, "SN_" + i_str + "_velocity", "0,0,0"));


  region_name = pin->GetOrAddString(block_name, "SN_" + i_str + "_region", "ball"); // region 的种类名，默认是 ball
  // 由于 radius 和 center 是比较通用的参数，因此总是读取
  radius = pin->GetOrAddReal(block_name, "SN_" + i_str + "_radius", 0.0);
  center_str = pin->GetOrAddString(block_name, "SN_" + i_str + "_center", "0,0,0"); // 用于后续的 fallback
  center = read_array(center_str);

  if (region_name == "ball") {
    // energy 的 region 是默认的，input file 里的 radius 和 center 都是指 energy 的。mass 的相关参数可以单独指定。
    energy_radius = radius;
    energy_center = center;
    energy_region = std::unique_ptr<Region>(new Ball(energy_center, energy_radius, dim));
    // mass 的 region。如果没有单独指定，则使用 SN_i_radius 和 SN_i_center 所指定的。
    // 目前不检查 mass 是否为 0，避免麻烦。内存开销应该很小。
    mass_radius = pin->GetOrAddReal(block_name, "SN_" + i_str + "_mass_radius", radius);
    mass_center = read_array(pin->GetOrAddString(block_name, "SN_" + i_str + "_mass_center", center_str));
    mass_region = std::unique_ptr<Region>(new Ball(mass_center, mass_radius, dim));
    
  } else if (region_name == "heart") {
    //TODO 把 heart SN 实现
    throw std::invalid_argument("Heart region is not implemented yet.");
  } else {
    throw std::invalid_argument("Unknown region type: " + region_name);
  }

  energy_density = energy / energy_region->volume; //* volume 有时会被初始化为 0，因为这种时候 region.contains() 永远不会为 true
  mass_density = mass / mass_region->volume;

  // ? 有没有必要检查 energy_radius <= mass_radius ?

}


std::vector<SuperNova> read_supernova_list(ParameterInput *pin, const int dim) {
  std::vector<SuperNova> supernova_list;
  std::string block_name, i_str;
  block_name = "supernova";
  for (int i = 1; true; i++) {
    i_str = std::to_string(i);
    if (pin->DoesParameterExist("supernova", "SN_" + i_str + "_time")) { // 目前根据是否有 SN_i_time 来判断是否有这个 SN
      supernova_list.push_back(SuperNova(pin, i, dim));
    } else {
      break;
    }
  }
  return supernova_list;
}


