#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

#include <cuda.h>

namespace stencils {
struct VelocityStencil {
  __device__ void operator()(FlowFieldView &field, int i, int j)
  {
    int obstacle = field.obstacles.get(i, j);
    // if fluid cell
    if ((obstacle & params::OBSTACLE_SELF) == 0) {
      if ((obstacle & params::OBSTACLE_RIGHT) == 0) {
        // if neighbour also fluid
        const auto dx = c_params.mesh.mesh_dx;
        field.v.u(i, j) =
            field.fgh.u(i, j) -
            d_timestep.dt / dx * (field.p.get(i + 1, j) - field.p.get(i, j));
      }
      else {
        // else set to 0
        field.v.u(i, j) = 0.;
      }
      if ((obstacle & params::OBSTACLE_TOP) == 0) {
        const auto dy = c_params.mesh.mesh_dy;
        field.v.v(i, j) =
            field.fgh.v(i, j) -
            d_timestep.dt / dy * (field.p.get(i, j + 1) - field.p.get(i, j));
      }
      else {
        field.v.v(i, j) = 0.;
      }
    }
  }
};
} // namespace stencils
