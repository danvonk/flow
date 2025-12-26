#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

#include <cuda.h>

namespace stencils {
struct ObstacleStencil {
  __device__ void operator()(FlowFieldView field, int i, int j)
  {
    const int obstacle = field.obstacles.get(i, j);
    const auto dx = c_params.mesh.mesh_dx;
    const auto dy = c_params.mesh.mesh_dy;
    // check if current cell is an obstacle
    if ((obstacle & params::OBSTACLE_SELF) == 1) {
      // if cell above is fluid then must enforce no-slip
      if ((obstacle & params::OBSTACLE_TOP) == 0) {
        field.v.u(i, j) = -dy / dy * field.v.u(i, j + 1);
      }
      if ((obstacle & params::OBSTACLE_BOTTOM) == 0) {
        field.v.u(i, j) = -dy / dy * field.v.u(i, j - 1);
      }
      if ((obstacle & params::OBSTACLE_RIGHT) == 0) {
        field.v.v(i, j) = -dx / dx * field.v.v(i + 1, j);
      }
      if ((obstacle & params::OBSTACLE_LEFT) == 0) {
        field.v.v(i, j) = -dx / dx * field.v.v(i - 1, j);
      }

      // set normal velocity to zero if right/top neighbour not obstacle
      if ((obstacle & params::OBSTACLE_RIGHT) == 0) {
        field.v.u(i, j) = 0.;
      }
      if ((obstacle & params::OBSTACLE_TOP) == 0) {
        field.v.v(i, j) = 0.;
      }
    }
  }
};

} // namespace stencils
