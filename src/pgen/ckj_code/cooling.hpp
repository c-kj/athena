#ifndef COOLING_HPP_
#define COOLING_HPP_

// C++ headers
// #include <cmath>

// Athena++ headers
// #include "../../athena.hpp"

void CoolingSourceTerm(std::string cooling_model, 
            MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar);





#endif // COOLING_HPP_