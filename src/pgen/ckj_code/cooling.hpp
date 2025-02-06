#ifndef COOLING_HPP_
#define COOLING_HPP_

// C++ headers
#include <memory>

// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"

// 自定义的头文件
#include "ckj_code.hpp"
#include "Townsend.hpp"

struct TownsendCooling; // forward declaration

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

struct Draine_2011_cutoff : public CoolingModel {
  Real CoolingCurve(Real T_cgs) const override;
};

struct OtherModel : public CoolingModel {
  Real CoolingCurve(Real T_cgs) const override;
  Real Jacobian(Real T_cgs) const override;
};


struct Cooling {
  Cooling() = default;  // 默认构造函数。需要这个才能在声明 Cooling cooling 时初始化。不过如果改用指针，可能就不需要这个了。
  Cooling(Mesh *pmy_mesh, ParameterInput *pin);

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
  bool use_prim_in_cooling;  // 使用 prim 或 cons 来计算 rho、P。如果用 prim，相当于一阶 Euler 积分；如果用 cons，相当于一阶 operator-splitting

  // Integrator
  std::string integrator;
  bool implicit;
  Real CFL_cooling;

  // Townsend cooling 机制的对象
  std::unique_ptr<TownsendCooling> Townsend;

  // RootFinder
  int max_iter;
  Real rel_tol, abs_tol;

  // Subcycle
  Real CFL_subcycle;
  bool subcycle_adaptive;

  // Limiter
  bool limiter_on;
  Real T_floor_cgs;

  // 检查 EOS 是 adiabatic 的，这样才有良定义的 gamma
  static_assert(NON_BAROTROPIC_EOS == 1, "### ERROR: Cooling 类目前只支持 NON_BAROTROPIC_EOS！因为要用到内能。");
  static_assert(GENERAL_EOS == 0, "### ERROR: Cooling 类目前不支持 GENERAL_EOS！因为要用到 gamma (=5/3)。");  // 若要使用 GENERAL_EOS，所有调用 gamma 的地方都需要修改！
  Real gamma;
  //FUTURE 目前这里不能声明为 const，否则需要在构造函数中初始化，就没法直接拷贝构造了。等把 cooling 改为用 unique_ptr 之后，可以改为 const。

  // 关于各种类粒子的丰度，所需要的参数，从 Abundance 类中获取。目前设为 static constexpr。
  static constexpr Real mu = Abundance::mu;
  static constexpr Real mu_e = Abundance::mu_e;
  static constexpr Real X_H = Abundance::X_H;

  
};



#endif // COOLING_HPP_