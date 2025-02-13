
#ifndef CKJ_CODE_REFINEMENT_CONDITION_HPP_
#define CKJ_CODE_REFINEMENT_CONDITION_HPP_

#include "../../athena.hpp"


// 自定义 AMR 机制所使用的变量
// 注意普通的变量是不能作为全局变量跨源文件的
   //TODO 目前用 extern 的方式定义全局变量，以后可以考虑改用类。
extern Real time_start_AMR;
extern std::vector<std::array<Real, 3>> point_list;
extern std::vector<int> level_list;

extern int AMR_grad_p, AMR_grad_rho;
extern Real grad_p_refine, grad_p_derefine, grad_rho_refine, grad_rho_derefine;

extern int my_root_level;     // Mesh 的 root_level 是 private 的，所以需要自己定义一个变量来存储它



int RefinementCondition_Point(MeshBlock *pmb, std::array<Real, 3> point, int level_limit);

int RefinementCondition_Gradient(MeshBlock *pmb);

int RefinementCondition(MeshBlock *pmb);

void InitRefinementCondition(ParameterInput *pin);

void SetRootLevel(int root_level);

#endif // CKJ_CODE_REFINEMENT_CONDITION_HPP_