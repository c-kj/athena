#ifndef TURB_BH_ACCRETION_HPP_
#define TURB_BH_ACCRETION_HPP_


// 自定义函数所需的标准库导入
#include <fstream>  // std::ofstream
#include <memory>   //std::unique_ptr

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

// 自定义的头文件
#include "ckj_code/supernova.hpp"
#include "ckj_code/cooling.hpp"
#include "ckj_code/initial_condition.hpp"


//BUG 这些全局变量并没有使用 extern 来声明，如果被多个源文件包含的话，会重复定义，违反 ODR，编译也会报错。目前这些全局变量只在 turb+BH_accretion.cpp 中使用，暂时没出问题。
//TODO 将这些全局变量大多数改造为类的成员变量。剩余的用 extern 声明，然后在 turb+BH_accretion.cpp 中定义。
//* 实际上，目前没有任何其他源文件（除了 turb+BH_accretion.cpp）包含这个 hpp。这里面的大部分声明都可以直接挪到 turb+BH_accretion.cpp 的开头。
// 声明自定义的全局变量，从 input file 中读取
// 这些变量的赋值是在 Mesh::InitUserMeshData 中完成的

Real GM_BH, R_in, R_out, rho_in_BH;
Real M_BH;




Supernovae supernovae;


Cooling cooling;  //* 暂时不用指针，而是直接（在栈上）创建一个对象。这样做的好处是，不用担心对象的生命周期问题。缺点是：对象的大小必须在编译时知道，不能动态改变。如果对象很大，可能会消耗大量的栈空间。

std::unique_ptr<InitialCondition> initial_condition;


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