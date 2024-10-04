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

Real GM_BH, R_in, R_out, dens_in_BH;
Real M_BH;


//FUTURE 这几个「自定义机制」对象/指针，考虑挪到 ckj_code.hpp 中去，从而可以在其他模块、pgen 中引用

Supernovae supernovae;

Cooling cooling;  //* 暂时不用指针，而是直接（在栈上）创建一个对象。这样做的好处是，不用担心对象的生命周期问题。缺点是：对象的大小必须在编译时知道，不能动态改变。如果对象很大，可能会消耗大量的栈空间。

std::unique_ptr<InitialCondition> initial_condition;


//TODO 把 debug 机制相关的变量也都挪到 utils 或 debugging 中去
int verbose, debug;
std::string debug_filepath, verbose_filepath;
std::ofstream debug_stream;



/* -------------------------------------------------------------------------- */
/*                                Source Terms                                */
/* -------------------------------------------------------------------------- */

// 本 pgen 中要 Enroll 的显式源项函数。在其中调用以下各个单独的源项函数。
//* 应该只依赖于 prim 而修改 cons。不能修改 prim。不应该依赖于 cons（尤其是因为 cons 可能被前面的源项修改）。
void MySourceFunction(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

// 只在 time integrator 的最后一个 stage 调用的源项函数。
void SourceTermAtLastStage(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

// 中心黑洞的引力
void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

// 中心黑洞作为 sink region 的处理：记录吸积相关物理量，并清空 sink 内的各种 cons。
//* 这个函数应该在所有的源项之后调用，因为吸积的量依赖于 cons，需要等其他源项对 cons 做完修改。
void SMBH_sink(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);





#endif // TURB_BH_ACCRETION_HPP_