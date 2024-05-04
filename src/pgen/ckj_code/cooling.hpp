#ifndef COOLING_HPP_
#define COOLING_HPP_

// C++ headers
// #include <cmath>

// Athena++ headers
// #include "../../athena.hpp"

//TODO
struct Cooling {
  Cooling() = default;  // 默认构造函数。需要这个才能在声明 Cooling cooling 时初始化。不过如果改用指针，可能就不需要这个了。
  Cooling(ParameterInput *pin, Mesh *pmy_mesh);
  
  const Real CoolingFunction(Real T_cgs);
  // Real Integrator();
  // void Limiter();

  const Real CoolingRate(const Real &rho, const Real &P);

  const Real CoolingTimeScale(const Real &rho, const Real &P);
  const Real CoolingTimeStep(MeshBlock *pmb);

  void CoolingSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                         const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                         const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                         AthenaArray<Real> &cons_scalar);

  Units *punit;
  std::string cooling_model, integrator;
  bool cooling_flag;
  bool operator_splitting, implicit;


  // Real gamma; //* 不用 gamma，以兼容 General EoS
  Real X_H, X_He, X_Metal;
  Real mu, mu_e;

  Real CFL_cooling;
  
};



#endif // COOLING_HPP_