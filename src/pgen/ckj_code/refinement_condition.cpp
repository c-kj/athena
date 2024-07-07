


// Athena++ headers
#include "../../athena.hpp"
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../../parameter_input.hpp"

// 自定义的头文件
#include "refinement_condition.hpp"
#include "utils.hpp"

// 自定义 AMR 机制所使用的变量

Real time_start_AMR;
std::vector<std::array<Real, 3>> point_list;
std::vector<int> level_list;

int AMR_grad_p, AMR_grad_rho;
Real grad_p_refine, grad_p_derefine, grad_rho_refine, grad_rho_derefine;

int my_root_level;     // Mesh 的 root_level 是 private 的，所以需要自己定义一个变量来存储它




// 通过判断 MeshBlock 是否包含指定的点，确定是否需要细化
int RefinementCondition_Point(MeshBlock *pmb, std::array<Real, 3> point, int level_limit=100000) {
  RegionSize size = pmb->block_size;
  int current_level = pmb->loc.level - my_root_level;
  if (current_level >= level_limit) return 0; // 如果当前 level 大于等于 level_limit，则不再细化

  Real x1 = point[0], x2 = point[1], x3 = point[2];
  if (size.x1min <= x1 && x1 <= size.x1max &&
      size.x2min <= x2 && x2 <= size.x2max &&
      size.x3min <= x3 && x3 <= size.x3max) {
    return 1;
  }
  return 0;
}


// TEMP 用于shock 的压强、密度梯度细化条件
int RefinementCondition_Gradient(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real dlnp = 0;
  Real dlnrho = 0;
  Real dlnp_new, dlnrho_new;
  if (AMR_grad_p) {
    if (pmb->pmy_mesh->f3) {
      for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
        for (int j=pmb->js-1; j<=pmb->je+1; j++) {
          for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
            dlnp_new = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                                +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                                +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
            dlnp = std::max(dlnp, dlnp_new);
          }
        }
      }
    } else if (pmb->pmy_mesh->f2) {
      int k = pmb->ks;
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          dlnp_new = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                              +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
          dlnp = std::max(dlnp, dlnp_new);
        }
      }
    }
  }

  if (AMR_grad_rho) {
    if (pmb->pmy_mesh->f3) {
      for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
        for (int j=pmb->js-1; j<=pmb->je+1; j++) {
          for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
            dlnrho_new = std::sqrt(SQR(0.5*(w(IDN,k,j,i+1) - w(IDN,k,j,i-1)))
                                  +SQR(0.5*(w(IDN,k,j+1,i) - w(IDN,k,j-1,i)))
                                  +SQR(0.5*(w(IDN,k+1,j,i) - w(IDN,k-1,j,i)))/w(IDN,k,j,i));
            dlnrho = std::max(dlnrho, dlnrho_new);
          }
        }
      }
    } else if (pmb->pmy_mesh->f2) {
      int k = pmb->ks;
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          dlnrho_new = std::sqrt(SQR(0.5*(w(IDN,k,j,i+1) - w(IDN,k,j,i-1)))
                                +SQR(0.5*(w(IDN,k,j+1,i) - w(IDN,k,j-1,i))))/w(IDN,k,j,i);
          dlnrho = std::max(dlnrho, dlnrho_new);
        }
      }
    }
  }
  // TEMP 临时注释掉 debug 机制，因为 debug 机制还没迁移到独立文件中
  // if (debug >= DEBUG_MeshBlock) {
  //   printf_and_save_to_stream(debug_stream, 
  //   "DEBUG_MeshBlock: dlnp = %f, dlnrho = %f, grad_p_refine = %f, grad_rho_refine = %f, AMR_grad_p = %d, AMR_grad_rho = %d \n", 
  //   dlnp, dlnrho, grad_p_refine, grad_rho_refine, AMR_grad_p, AMR_grad_rho);
  // }

  if (dlnp > grad_p_refine || dlnrho > grad_rho_refine) return 1;
  if (dlnp < grad_p_derefine && dlnrho < grad_rho_derefine) return -1;
  return 0;
}


//TODO 把这个函数改名，避免跟其他 pgen 中的撞了？
//? 是否有必要把函数的声明和定义分开？声明放在前面，定义放在后面？
// 测试能不能用 AMR 来充当 SMR，初始化湍流
int RefinementCondition(MeshBlock *pmb) {
  // 在 time_start_AMR 之前，都不进行 AMR，直接返回 0
  Real time = pmb->pmy_mesh->time;
  if (time < time_start_AMR) return 0;

  // TEMP
  // debug 消息，放在前面，免得还没输出就 return 了
  // if (debug >= DEBUG_Mesh && pmb->gid == 0) {
  //   printf_and_save_to_stream(debug_stream, "DEBUG_Mesh: RefinementCondition is called at time = %f \n", time);
  // }

  //* 这里涉及有些复杂的逻辑判断：每一个 condition 返回的是 -1,0,1 之一，如果有 -1 和 1 的冲突那么需要仔细考量优先级。
  for (int i = 0; i < point_list.size(); i++) {
    if (RefinementCondition_Point(pmb, point_list[i], level_list[i]) == 1) return 1; // 如果 MeshBlock 包含任何一个指定的点，则进行细化
  }
  // 这里这么写是为了后续扩展
  // TEMP
  return RefinementCondition_Gradient(pmb);

  return 0;
}

//TODO 把这里的 INFINITY 改为 std::numeric_limits<Real>::max()
void InitRefinementCondition(ParameterInput *pin) {
  // 读取自定义的 AMR 参数
  time_start_AMR = pin->GetOrAddReal("AMR", "time_start_AMR", 0.0);

  std::tie(point_list, level_list) = get_AMR_points_and_levels(pin);  // 用 tie 函数将返回的 tuple 解包

  AMR_grad_p = pin->DoesParameterExist("AMR/gradient", "grad_p_refine") || pin->DoesParameterExist("AMR/gradient", "grad_p_derefine");
  AMR_grad_rho = pin->DoesParameterExist("AMR/gradient", "grad_rho_refine") || pin->DoesParameterExist("AMR/gradient", "grad_rho_derefine");
  grad_p_refine = pin->GetOrAddReal("AMR/gradient", "grad_p_refine", INFINITY);
  grad_p_derefine = pin->GetOrAddReal("AMR/gradient", "grad_p_derefine", -INFINITY);
  grad_rho_refine = pin->GetOrAddReal("AMR/gradient", "grad_rho_refine", INFINITY);
  grad_rho_derefine = pin->GetOrAddReal("AMR/gradient", "grad_rho_derefine", -INFINITY);
  // TEMP
  // if (debug >= DEBUG_Main) {
  //   printf_and_save_to_stream(debug_stream, "DEBUG_Main: AMR_grad_p = %d, AMR_grad_rho = %d \n", AMR_grad_p, AMR_grad_rho);
  // }

}

void SetRootLevel(int root_level) {
  my_root_level = root_level;
}


