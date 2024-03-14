#ifndef TURB_BH_ACCRETION_HPP_
#define TURB_BH_ACCRETION_HPP_


// 自定义函数所需的标准库导入
#include <fstream>  // std::ofstream

// Athena++ headers
// #include "../athena.hpp"
// #include "../athena_arrays.hpp"
// #include "../coordinates/coordinates.hpp"
// #include "../eos/eos.hpp"
// #include "../fft/athena_fft.hpp"
// #include "../field/field.hpp"
// #include "../hydro/hydro.hpp"
// #include "../parameter_input.hpp"
// #include "../utils/utils.hpp"

#include "ckj_code/supernova.hpp"


// 声明自定义的全局变量，从 input file 中读取
// 这些变量的赋值是在 Mesh::InitUserMeshData 中完成的
// 使用 namespace 创建一个命名空间 my_vars，然后再 using namespace 导入它
// 这样做的好处是，所有自定义的变量都放在一起，在大纲层级中更加清晰

Real GM_BH, R_in, R_out, rho_in_BH, rho_init, E_tot_init;





std::string init_cond_type;
Real rho_at_boundary, E_tot_at_boundary, power_law_index, alpha;
Real approx_Bondi_rho_profile(Real alpha, Real R_Bondi, Real r);

//TEMP
int SN_flag;
std::string SN_type;
Real SN_time, SN_energy, SN_radius, SN_mass;

std::vector<SuperNova> supernova_list;
std::vector<SuperNova*> supernova_to_inject;



//TODO 把 debug 机制相关的变量也都挪到 utils 或 debugging 中去
int verbose, debug;
std::string debug_filepath, verbose_filepath;
std::ofstream debug_stream;




// 源项
void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);





#endif // TURB_BH_ACCRETION_HPP_