#include <limits>

#include "initial_condition.hpp"



// 构造函数
InitialCondition::InitialCondition(Mesh *pmesh, ParameterInput *pin): 
    pin{pin},
    punit{pmesh->punit},
    gamma{pin->GetReal("hydro", "gamma")}, 
    init_cond_type{pin->GetOrAddString("initial_condition","init_cond_type","uniform")},
    n_init_cgs{pin->GetReal("initial_condition", "n_init_cgs")},
    T_init_cgs{pin->GetReal("initial_condition", "T_init_cgs")}
{
  // 将 n 和 T 换算为 code unit 下的 rho 和 E_thermal。在后续使用时，以这两个量作为基础。
  rho_init_code = rho_from_n_cgs(mu, n_init_cgs) / punit->code_density_cgs;
  E_thermal_init_code = P_from_nT_cgs(n_init_cgs, T_init_cgs) / (gamma - 1.0) / punit->code_energydensity_cgs;
}

// void InitialCondition::PrintInfo() {
//   std::cout << "Initial condition type: " << init_cond_type << "\n";
//   std::cout << "n_init_cgs = " << n_init_cgs << "\n";
//   std::cout << "T_init_cgs = " << T_init_cgs << "\n";
//   std::cout << "rho_init_code = " << rho_init_code << "\n";
//   std::cout << "E_thermal_init_code = " << E_thermal_init_code  << "\n";
// }

// 遍历每个 cell，调用 SetSingleCell 来设定初始条件
void InitialCondition::SetInitialCondition(MeshBlock *pmb) {
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        SetSingleCell(pmb, i, j, k);
      }
    }
  }
};

// 对单个 cell 设定初始条件
inline void InitialCondition::SetSingleCell(MeshBlock *pmb, const int i, const int j, const int k) {
  auto & cons = pmb->phydro->u;
  // 计算坐标，因为大多数初始条件都是坐标的函数
  const Real x3 = pmb->pcoord->x3v(k);
  const Real x2 = pmb->pcoord->x2v(j);
  const Real x1 = pmb->pcoord->x1v(i);

  Real rho, vx1, vx2, vx3, E_thermal; //* 这五个量需要在 init_cond_type 的各个分支中赋值，最后转换成 cons

  // 均匀初始条件（默认）
  if (init_cond_type == "uniform") {
    rho = rho_init_code;
    //TODO 均匀初始速度，由速度的大小和方向（未归一化的矢量）设定
    vx1 = 0.0;
    vx2 = 0.0;
    vx3 = 0.0;
    E_thermal = E_thermal_init_code;
  }

  //TODO 其他的初始条件有待测试
  // 幂律初值，密度和能量遵循相同的幂律，初速度为 0。把 init 的值理解为 R_out 处的值。
  else if (init_cond_type == "power_law") {
    static Real R_out = pin->GetReal("initial_condition", "R_out");  //? 从 initial_condition 这个 block 获取吗？
    static Real power_law_index = pin->GetReal("initial_condition", "power_law_index"); // 对于类 Bondi profile，一般为负值

    Real r = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
    Real factor = pow(r/R_out, power_law_index);
    rho = rho_init_code * factor; 
    //? power law 的 径向速度？
    vx1 = 0.0;
    vx2 = 0.0;
    vx3 = 0.0;
    E_thermal = E_thermal_init_code * factor;
  }

  // approximate_Bondi: 初速度 = 0，密度和能量遵循 Bondi profile 的近似解
  else if (init_cond_type == "approximate_Bondi") {
    static Real velocity_factor = pin->GetOrAddReal("initial_condition","velocity_factor", 1.0);
    static Real alpha = pin->GetOrAddReal("initial_condition","alpha",1./3.);
    static Real R_out = pin->GetOrAddReal("problem", "R_out", std::numeric_limits<Real>::max());
    static Real M_BH = pin->GetOrAddReal("problem", "M_BH", 0.0) * punit->solar_mass_code;
    static Real GM_BH = M_BH * punit->grav_const_code;
    static Real c_s = sqrt(gamma * (gamma - 1.0) * E_thermal_init_code / rho_init_code );  //? 这里算声速，应该使用 gamma 还是 polytropic_index?
    static Real R_Bondi = 2 * GM_BH / SQR(c_s);
    // 目前这里是把 init 的值直接理解为边界值，然后换算出无穷远的值。以后可以考虑修改？
    static Real rho_inf = rho_init_code / approx_Bondi_rho_profile(alpha, R_Bondi, R_out);

    static Real polytropic_index = pin->GetOrAddReal("initial_condition","polytropic_index", gamma);  // 如果 input 文件中设置了 polytropic_index，那么就用它，否则用 gamma

    // 用于计算 Bondi 吸积率中的因子的函数。这里自变量 gamma 实际上是 polytropic_index
    auto M_dot_factor = [](Real gamma) -> Real {
      if (gamma >= 1.6666) return 1.0;           // gamma 接近 5/3 的情况下，需要取极限，M_dot_factor = 1
      if (gamma == 1.0   ) return exp(3./2);  // gamma == 1 的情况下，需要取极限，M_dot_factor = exp(3/2)
      return pow(2/(5 - 3*gamma),(5 - 3*gamma)/(2.*(gamma - 1)));
    };
    //BUG 这里使用的是 hydro 的 gamma，但若开启 cooling，实际上可能 isothermal 更合适。
    //TODO 所以，考虑在 cooling_on 的情况下使用 block 内自定义的 polytropic_index 参数。但要弄清楚 c_s 的计算中使用的 gamma 应该用哪个。
    static Real M_dot = PI * rho_inf * SQR(GM_BH) / CUBE(c_s) * M_dot_factor(polytropic_index);  

    Real r = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
    rho = rho_inf * approx_Bondi_rho_profile(alpha, R_Bondi, r);
    //? approximate Bondi 的 径向速度？
    Real v = M_dot / (4 * PI * SQR(r) * rho);  // 根据连续性方程计算出速度
    v *= velocity_factor;  // 额外乘一个系数
    vx1 = -v * x1 / r;
    vx2 = -v * x2 / r;
    vx3 = -v * x3 / r;
    E_thermal = E_thermal_init_code * pow(rho/rho_init_code, polytropic_index);
  } 

  // heart: 在心形区域内，n 和 T 设为与外部不同的值
  else if (init_cond_type == "heart") {
    std::string block = "initial_condition/heart";
    static Real a = pin->GetOrAddReal(block, "a", 3.3);
    static Real b = pin->GetOrAddReal(block,"b", 0.75);
    static Real heart_size = pin->GetOrAddReal(block,"heart_size", 1.0);

    static Real n_inside = pin->GetReal(block,"n_inside");
    static Real T_inside = pin->GetReal(block,"T_inside");

    static Real rho_inside = rho_from_n_cgs(mu, n_inside) / punit->code_density_cgs;
    static Real E_thermal_inside = P_from_nT_cgs(n_inside, T_inside) / (gamma - 1.0) / punit->code_energydensity_cgs;

    static Heart heart({0, 0, 0}, heart_size, a, b); // 构造一个心形区域。目前默认 center 是 (0,0,0)

    if (heart.contains({x1,x2,x3})) {
      rho = rho_inside;
      E_thermal = E_thermal_inside;
    } else {  // 在心形区域外
      rho = rho_init_code;
      E_thermal = E_thermal_init_code;
    }
    vx1 = 0.0;
    vx2 = 0.0;
    vx3 = 0.0;
  // } else if () {  // 其他的初始条件
  } else {
    throw std::runtime_error("Unknown init_cond_type !");
  }

  // 如果前面哪里忘了给 rho, v, E_thermal 赋值，那么这里 linter 就会提示说有 garbage value （但编译不会给 warning）
  cons(IDN,k,j,i) = rho;
  cons(IM1,k,j,i) = rho * vx1;
  cons(IM2,k,j,i) = rho * vx2;
  cons(IM3,k,j,i) = rho * vx3;
  if (NON_BAROTROPIC_EOS) {
    cons(IEN,k,j,i) = E_thermal + 0.5 * rho * (SQR(vx1) + SQR(vx2) + SQR(vx3));
  }

  // 如果启用 Passive Scalars，那么需要初始化它们。
  //TODO 目前全部初始化为 0。将来可能根据特定的 Passive Scalar 有不同的初始化。
  if (NSCALARS > 0) {
    for (int n=0; n<NSCALARS; ++n) {
      pmb->pscalars->s(n,k,j,i) = 0;
    }
  }
}

