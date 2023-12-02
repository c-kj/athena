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

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void SMBH_grav(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  // 从 MeshBlock 中读取参数
  Real GM = pmb->ruser_meshblock_data[0](0); 
  Real R_in = pmb->ruser_meshblock_data[0](1);
  Real R_out = pmb->ruser_meshblock_data[0](2);
  Real rho_in_BH = pmb->ruser_meshblock_data[0](3);
  Real rho_init = pmb->ruser_meshblock_data[0](4);
  Real E_tot_init = pmb->ruser_meshblock_data[0](5);

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
          cons(IM1,k,j,i) += - GM*x/r3 * rho * dt;
          cons(IM2,k,j,i) += - GM*y/r3 * rho * dt;
          cons(IM3,k,j,i) += - GM*z/r3 * rho * dt;

          // 能量的改变
          // ? 我看 Athena++ 的源代码都是判断了 if NON_BAROTROPIC_EOS，但这玩意默认是 True，不知道啥意思… 绝热方程是不是 barotropic？从模拟结果上来看好像不是
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) += - GM*(x*vx + y*vy + z*vz)/r3 * rho * dt; // 这里可以稍微优化，快一点
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


        if ( std::isnan(cons(IDN,k,j,i)) ) {
          printf("rho is nan at x,y,z,r = (%f, %f, %f, %f) \n", x, y, z, r);
        }

      }
    }
  }
  // std::cout << "my print: time = " << time << std::endl; // 从 print 的结果看，似乎每个时间步会调用两次，一次在开始，一次在中间
  return;
}



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

  // TODO 
EnrollUserExplicitSourceFunction(SMBH_grav);

  return;
}


// MODIFIED 读取 input 参数
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // 从 input file 中读取，放到 MeshBlock 中的数组中，这样就可以传递给自定义的源项了
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(10);
  ruser_meshblock_data[0](0) = pin->GetReal("problem","GM_BH");
  ruser_meshblock_data[0](1) = pin->GetReal("problem","R_in");
  ruser_meshblock_data[0](2) = pin->GetReal("problem","R_out");
  ruser_meshblock_data[0](3) = pin->GetReal("problem","rho_in_BH");
  ruser_meshblock_data[0](4) = pin->GetReal("problem","rho_init");
  ruser_meshblock_data[0](5) = pin->GetReal("problem", "E_tot_init");
  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = 1.0;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          Real E_tot_init = pin->GetReal("problem", "E_tot_init");
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
}
