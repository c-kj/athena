
#include <sstream>  // std::stringstream, std::getline
#include <stdexcept>  // std::invalid_argument

// #include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"


#include "supernova.hpp"
#include "utils.hpp"


// SuperNova 的构造函数
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

  // 这里 energy_density 和 mass_density 的单位有点奇怪：energy 是 input file 中的数字，单位为 1e51 erg（目前），而 volume 是 code unit。但这样做是为了尽量避免重复计算 / 重复储存。
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


// SuperNovae 的构造函数
SuperNovae::SuperNovae(Mesh *pmy_mesh, ParameterInput *pin) : punit(pmy_mesh->punit), ndim(pmy_mesh->ndim) {
  SN_flag = pin->GetOrAddInteger("supernova","SN_flag",0);  //* 目前 SN_flag 设为一个 int，是为了以后可能会有多种 SN 的情况（暂时也没想出来）
  if (SN_flag > 0) {
    if ( !NON_BAROTROPIC_EOS ) {
      throw std::runtime_error("SN explosion is not implemented for barotropic EOS!"); 
    }

    // 这里暂时没有写开关：什么情况下把 SN 的 input 单位解读为 1e51 erg 和 Msun。若有需求再改
    SN_energy_unit = punit->bethe_code;
    SN_mass_unit = punit->solar_mass_code;

    integrator = pin->GetOrAddString("time", "integrator", "vl2"); // 记录 hydro 所使用的 integrator。目前用来处理超新星的注入时机
    
    supernova_list = read_supernova_list(pin, ndim); // 需要传入 ndim，从而确定区域的维数
  }
}


void SuperNovae::GetSupernovaeToInject(const Real time, const Real dt) {
  // 找到当前时间步内需要注入的超新星，放入 supernova_to_inject 列表
  supernova_to_inject = {}; // 每次都要先清空
  for (SuperNova& SN : supernova_list) {
    for (Real SN_time : SN.time_list) {
      if (time <= SN_time && SN_time < time + dt) {
        supernova_to_inject.push_back(&SN); // &SN 解引用，取 SN 的地址，也就是一个指向 SN 的指针
      }
    }
  }
}

void SuperNovae::SuperNovaeSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar) {

  const Real mesh_time = pmb->pmy_mesh->time;
  const Real mesh_dt = pmb->pmy_mesh->dt;

  // 当前调用 SourceTerm 时的注入比例。
  // 对于 AfterSourceTerm 和 UserWorkInLoop，传入的 dt 都是 mesh_dt，从而 inject_ratio = 1.0
  // 对于 InSourceTerm，传入的 dt 是当前 stage 的 dt，inject_ratio = beta_last_stage，根据当前 stage 的比例来调整注入的比例。
  // 注意：传入到这个函数的 dt 不一定是 总 SourceTerm 函数在当前 stage 接收到的 dt。因此这个变量不叫 dt_ratio，以示区分。
  const Real inject_ratio = dt/mesh_dt;

  GetSupernovaeToInject(mesh_time, mesh_dt); // 寻找应在当前步中爆炸的 SN 时，使用 mesh_time 和 mesh_dt，而非 SourceTerm 传入的 time 和 dt
  if ( supernova_to_inject.empty() ) return; // 如果当前时间步没有超新星爆炸，直接跳出
  
  Real x,y,z,r; 
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
#pragma omp simd  // OpenMP 的 SIMD 并行。在纯 MPI 并行时可能不起作用？
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);

//TODO 限制不能在黑洞内部注入超新星！
        if (true) {
         // 超新星爆炸的能量和质量注入。//* 目前没有实现动量注入，也没有考虑 SN 的运动速度
          for (SuperNova* SN : supernova_to_inject) {
            if ( SN->energy_region->contains({x,y,z}) ) {
              cons(IEN,k,j,i) += SN->energy_density * SN_energy_unit * inject_ratio;
            }
            if ( SN->mass_region->contains({x,y,z}) ) {
              cons(IDN,k,j,i) += SN->mass_density * SN_mass_unit * inject_ratio;
            }
            // 这里的 debug 信息不再那么有用，将来可以考虑删去
            // if (debug >= DEBUG_Main && pmb->gid == 0 && SN_flag > 0){
            //   if(mesh_time <= SN_time && SN_time < mesh_time + mesh_dt && dt == mesh_dt ) {
            //     printf_and_save_to_stream(debug_stream, 
            //     "DEBUG_Mesh: SN explosion at time = %f, mesh_time = %f, dt = %f, mesh_dt = %f \n", 
            //     time, pmb->pmy_mesh->time, dt, pmb->pmy_mesh->dt);
            //   }
            // } 
          }
        }
      }
    }
  }
  return;
}
