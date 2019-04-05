//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file boundary_flag.cpp
//  \brief utilities for processing the user's input <mesh> ixn_bc, oxn_bc parameters and
// the associated internal BoundaryFlag enumerators

// C headers

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn GetBoundaryFlag(std::string input_string)
//  \brief Parses input string to return scoped enumerator flag specifying boundary
//  condition. Typically called in Mesh() ctor and in pgen/*.cpp files.

BoundaryFlag GetBoundaryFlag(const std::string& input_string) {
  if (input_string == "reflecting") {
    return BoundaryFlag::reflect;
  } else if (input_string == "outflow") {
    return BoundaryFlag::outflow;
  } else if (input_string == "user") {
    return BoundaryFlag::user;
  } else if (input_string == "periodic") {
    return BoundaryFlag::periodic;
  } else if (input_string == "shear_periodic") {
    return BoundaryFlag::shear_periodic;
  } else if (input_string == "polar") {
    return BoundaryFlag::polar;
  } else if (input_string == "polar_wedge") {
    return BoundaryFlag::polar_wedge;
  } else if (input_string == "none") {
    return BoundaryFlag::undef;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryFlag" << std::endl
        << "Input string=" << input_string << "\n"
        << "is an invalid boundary type" << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn GetBoundaryString(BoundaryFlag input_flag)
//  \brief Parses enumerated type BoundaryFlag internal integer representation to return
//  string describing the boundary condition. Typicall used to format descriptive errors
//  or diagnostics. Inverse of GetBoundaryFlag().

std::string GetBoundaryString(BoundaryFlag input_flag) {
  if (input_flag == BoundaryFlag::reflect) {
    return "reflecting";
  } else if (input_flag == BoundaryFlag::outflow) {
    return "outflow";
  } else if (input_flag == BoundaryFlag::user) {
    return "user";
  } else if (input_flag == BoundaryFlag::periodic) {
    return "periodic";
  } else if (input_flag == BoundaryFlag::shear_periodic) {
    return "shear_periodic";
  } else if (input_flag == BoundaryFlag::polar) {
    return "polar";
  } else if (input_flag == BoundaryFlag::polar_wedge) {
    return "polar_wedge";
  } else if (input_flag == BoundaryFlag::undef) {
    return "none";
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryString" << std::endl
        << "Input enum class BoundaryFlag=" << static_cast<int>(input_flag) << "\n"
        << "is an invalid boundary type" << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn CheckBoundaryFlag(BoundaryFlag block_flag, CoordinateDirection dir)
//  \brief Called in each MeshBlock's BoundaryValues() constructor. Mesh() ctor only
//  checks the validity of user's input mesh/ixn_bc, oxn_bc string values corresponding to
//  a BoundaryFlag enumerator before passing it to a MeshBlock and then BoundaryBase
//  object. However, not all BoundaryFlag enumerators can be used in all directions as a
//  valid MeshBlock boundary.

void CheckBoundaryFlag(BoundaryFlag block_flag, CoordinateDirection dir) {
  std::stringstream msg;
  msg << "### FATAL ERROR in CheckBoundaryFlag" << std::endl
      << "Attempting to set invalid MeshBlock boundary= " << GetBoundaryString(block_flag)
      << "\nin x" << dir+1 << " direction" << std::endl;
  switch(dir) {
    case CoordinateDirection::X1DIR:
      switch(block_flag) {
        case BoundaryFlag::polar:
        case BoundaryFlag::polar_wedge:
        case BoundaryFlag::undef:
          ATHENA_ERROR(msg);
          break;
        default:
          break;
      }
      break;
    case CoordinateDirection::X2DIR:
      switch(block_flag) {
        case BoundaryFlag::shear_periodic:
        case BoundaryFlag::undef:
          ATHENA_ERROR(msg);
          break;
        default:
          break;
      }
      break;
    case CoordinateDirection::X3DIR:
      switch(block_flag) {
        case BoundaryFlag::polar:
        case BoundaryFlag::polar_wedge:
        case BoundaryFlag::shear_periodic:
        case BoundaryFlag::undef:
          ATHENA_ERROR(msg);
          break;
        default:
          break;
      }
  }
  return;
}
