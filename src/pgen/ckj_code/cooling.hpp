#ifndef COOLING_HPP_
#define COOLING_HPP_

// C++ headers
#include <memory>

// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"

// 自定义的头文件
#include "ckj_code.hpp"

// CoolingModel 是各种冷却函数的基类
//* 把 CoolingModel 写成类，本来是为了实现对应的 Jacobian。但目前 Jacobian 并没有用到，所以也没啥意义……
class CoolingModel {
public:
  virtual ~CoolingModel() {}

  // 返回冷却函数的值 (in cgs unit)
  Real CoolingFunction(Real T_cgs) const {
    if (T_cgs <= 0.0) { return 0.0;} // 避免 T<0 时指数截断反而变得巨大
    return CoolingCurve(T_cgs);
  }

  virtual Real CoolingCurve(Real T_cgs) const = 0;  // CoolingCurve 不应直接使用，因为它没有对输入的 T 做检查。
  virtual Real Jacobian(Real T_cgs) const = 0;  //* 目前用不上这个 Jacobian，而是直接都从有限差分获得。

  static std::unique_ptr<CoolingModel> Create(const std::string &cooling_model);

protected:
  const Real min_dT = std::sqrt(std::numeric_limits<Real>::epsilon());
};

class Draine_2011 : public CoolingModel {
public:
  Real CoolingCurve(Real T_cgs) const override;

  Real Jacobian(Real T_cgs) const override;
};

class Model2 : public CoolingModel {
public:
  Real CoolingCurve(Real T_cgs) const override;

  Real Jacobian(Real T_cgs) const override;
};


struct Cooling {
  Cooling() = default;  // 默认构造函数。需要这个才能在声明 Cooling cooling 时初始化。不过如果改用指针，可能就不需要这个了。
  Cooling(ParameterInput *pin, Mesh *pmy_mesh);
  
  // Real CoolingFunction(Real T_cgs) const;

  Real Integrator(const std::function<Real(Real)>& RHS, Real y0, Real dt) const;
  static Real RootFinder(const std::function<Real(Real)>& func, Real y0, int max_iter, Real rtol, Real atol);
  static Real FiniteDifferenceDerivative(const std::function<Real(Real)>& func, Real x);
  // void Limiter();

  Real CoolingRate(const Real &rho, const Real &P) const;

  Real CoolingTimeScale(const Real &rho, const Real &P) const;
  Real CoolingTimeStep(MeshBlock *pmb) const;

  void CoolingSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                         const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                         const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                         AthenaArray<Real> &cons_scalar);
  Units *punit;

  // CoolingModel
  std::string cooling_model;
  std::unique_ptr<CoolingModel> model;
  bool cooling_flag;

  // Cooling SourceTerm 的注入方式
  SourceTermPosition source_term_position;
  bool use_prim_in_cooling;

  // Integrator
  std::string integrator;
  bool implicit;
  Real CFL_cooling;

  // RootFinder
  int max_iter;
  Real rel_tol, abs_tol;

  // Subcycle
  Real CFL_subcycle;
  bool subcycle_adaptive;

  // Limiter
  bool limiter_on;
  Real T_floor_cgs;

  // Real gamma; //* 不用 gamma，以兼容 General EoS
  Real X_H, X_He, X_Metal;
  Real mu, mu_e;

  
};



#endif // COOLING_HPP_