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
struct CoolingModel {
  virtual ~CoolingModel() = default;

  static std::unique_ptr<CoolingModel> Create(const std::string &cooling_model);

  // 返回冷却函数的值 (in cgs unit)
  Real CoolingFunction(Real T_cgs) const {
    if (T_cgs <= 0.0) { return 0.0;} // 避免 T<0 时指数截断反而变得巨大
    return CoolingCurve(T_cgs);
  }

  virtual Real CoolingCurve(Real T_cgs) const = 0;  // CoolingCurve 不应直接使用，因为它没有对输入的 T 做检查。

  // Jacobian 是虚函数但不是纯虚函数，这里提供默认实现
  virtual Real Jacobian(Real T_cgs) const {
    throw std::runtime_error("目前并不需要使用 CoolingModel 的 Jacobian！");  // 因为真正需要计算的导数是在 Newton-Raphson 求根时，在那里计算导数更方便，无需考虑复合函数的导数转换
#if false
    auto cooling_curve = [this](Real T_cgs) { return CoolingCurve(T_cgs); };
    return Cooling::FiniteDifferenceDerivative(cooling_curve, T_cgs);
#endif
  };
};

struct Draine_2011 : public CoolingModel {
  Real CoolingCurve(Real T_cgs) const override;
};

struct OtherModel : public CoolingModel {
  Real CoolingCurve(Real T_cgs) const override;
  Real Jacobian(Real T_cgs) const override;
};


struct Cooling {
  Cooling() = default;  // 默认构造函数。需要这个才能在声明 Cooling cooling 时初始化。不过如果改用指针，可能就不需要这个了。
  Cooling(ParameterInput *pin, Mesh *pmy_mesh);

  Real Integrator(const std::function<Real(Real)>& RHS, Real y0, Real dt) const;
  static Real RootFinder(const std::function<Real(Real)>& func, Real x0, int max_iter, Real rel_tol, Real abs_tol);
  static Real FiniteDifferenceDerivative(const std::function<Real(Real)>& func, Real x);
  // void Limiter();

  Real CoolingRate(const Real rho, const Real P) const;

  static Real CoolingTimeScale(const Real E_thermal, const Real cooling_rate);
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