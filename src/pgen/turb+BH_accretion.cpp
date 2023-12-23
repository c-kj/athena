//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//! \brief Problem generator for turbulence driver

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// 声明自定义的全局变量，从 input file 中读取
// 这些变量的赋值是在 Mesh::InitUserMeshData 中完成的
// 使用 namespace 创建一个匿名的命名空间。这会使得这些变量只能在本文件中使用，不会和其他文件中的同名变量冲突
// 这样做的好处是，所有自定义的变量都放在一起，在大纲层级中更加清晰
namespace {
  Real GM_BH, R_in, R_out, rho_in_BH, rho_init, E_tot_init;
  Real current_time, time_start_AMR;
  int verbose;

}


//? 是否有必要把函数的声明和定义分开？声明放在前面，定义放在后面？
// 测试能不能用 AMR 来充当 SMR，初始化湍流
int RefinementCondition(MeshBlock *pmb) {
  // 在 time_start_AMR 之前，都不进行 AMR，直接返回 0
  if (current_time < time_start_AMR) return 0;

  Real x,y,z,dx,dy,dz;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    dz = pmb->pcoord->dx3f(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      dy = pmb->pcoord->dx2f(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        dx = pmb->pcoord->dx1f(i);

        // 如果当前 MeshBlock 中有任何一个格点，其 xyz 坐标绝对值均小于其格子边长，则认定该 MeshBlock 紧贴坐标原点，需要细化
        if (std::abs(x) < dx && std::abs(y) < dy && std::abs(z) < dz) return 1;
      }
    }
  }
  
  return 0;
}



void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  current_time = time;

  Real x,y,z,r,r3,vx,vy,vz; 
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        r = sqrt(x*x + y*y + z*z);
        r3 = r*r*r;
        
        vx = prim(IVX,k,j,i);
        vy = prim(IVY,k,j,i);
        vz = prim(IVZ,k,j,i);
        Real& rho = cons(IDN,k,j,i); // rho 是引用，修改 rho 会修改 cons(IDN,k,j,i) 的值

        // TODO 只对黑洞外的区域施加引力，这样处理正确吗？
        // TODO 这里需要 && r <= R_out 吗？
        if (r >= R_in ) {
          cons(IM1,k,j,i) += - GM_BH*x/r3 * rho * dt;
          cons(IM2,k,j,i) += - GM_BH*y/r3 * rho * dt;
          cons(IM3,k,j,i) += - GM_BH*z/r3 * rho * dt;

          // 能量的改变
          // ? 我看 Athena++ 的源代码都是判断了 if NON_BAROTROPIC_EOS，但这玩意默认是 True，不知道啥意思… 绝热方程是不是 barotropic？从模拟结果上来看好像不是
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) += - GM_BH*(x*vx + y*vy + z*vz)/r3 * rho * dt; // 这里可以稍微优化，快一点
          }

        }

        // 进入黑洞的处理
        if (r < R_in) {
          // rho = std::min(rho_floor, rho); // 进入黑洞后，密度取 min(rho_floor, rho) // TODO 这个正确吗？有必要吗？可能自带约束？
          rho = rho_in_BH;

          // TODO 修改动量使其不外流。这样做正确吗？
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;

          // 修改能量
          cons(IEN,k,j,i) = 0.0 + 0.0; // TODO 应该怎么修改能量？尤其是内能怎么计算？ shocks?

        }

        // TODO 外部区域的处理，应该怎么做？
        // 目前是设定为保持初始值
        // TODO 或许可以改为继承该格子之前的值？但怎么操作呢？在 source term 里面不知道能不能做这件事
        if (r > R_out) {
          rho = rho_init;
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          cons(IEN,k,j,i) = E_tot_init;
        }

        // debug 信息：逐个格点检查 rho 是否为 nan
        if (verbose >= 5) {
          if ( std::isnan(cons(IDN,k,j,i)) ) {
            printf("rho is nan at x,y,z,r = (%f, %f, %f, %f) \n", x, y, z, r);
          }
        }

      }
    }
  }
    if (verbose >= 2 && pmb->gid == 0 ) {
      printf("source term at time = %f \n", time); // 从 print 的结果看，似乎每个时间步会调用两次，一次在开始，一次在中间
    }
  return;
}





//========================================================================================
// 以下是 Athena++ 提供的接口


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem","four_pi_G");
    SetFourPiG(four_pi_G);
  }
  
  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for impulsively driven turbulence
  // turb_flag = 3 for continuously driven turbulence
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  // 从 input file 中读取参数，存到 pgen 文件的全局变量中
  GM_BH = pin->GetReal("problem","GM_BH");
  R_in = pin->GetReal("problem","R_in");
  R_out = pin->GetReal("problem","R_out");
  rho_in_BH = pin->GetReal("problem","rho_in_BH");
  rho_init = pin->GetReal("problem","rho_init");
  E_tot_init = pin->GetReal("problem", "E_tot_init");

  time_start_AMR = pin->GetOrAddReal("problem", "time_start_AMR", 0.0);
  verbose = pin->GetOrAddInteger("problem", "verbose", 0);

  // 将自定义的源项注册到 Athena++ 中
  EnrollUserExplicitSourceFunction(SMBH_grav);

  // 将自定义的 AMR 条件注册到 Athena++ 中
  if(adaptive){
    EnrollUserRefinementCondition(RefinementCondition);
  }

  return;
}


// 这里可以设定属于每个 MeshBlock 的自定义数组数据。之前我用来传递 input 参数，后来改成了用 pgen 的全局变量
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

// MeshBlock::ProblemGenerator 用于设定初始条件（每个 MeshBlock 的）
// defined in either the prob file or default_pgen.cpp in ../pgen/
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho_init;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = E_tot_init;
        }
      }
    }
  }
}


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // 这个函数用于在模拟结束后做一些事情
  if (verbose >= 2) {
    printf("Mesh::UserWorkAfterLoop is called \n");
  }
}
