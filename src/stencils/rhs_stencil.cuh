#pragma once

#include "config.hpp"
#include "fields/flow_field.hpp"
#include "parameters.hpp"

#include <cuda.h>

namespace stencils {
class RHSStencil {
public:
  __device__ void operator()(FlowFieldView &field, int i, int j) const
  {
    field.rhs.get(i, j) =
        (1.0 / d_timestep.dt) *
        ((field.fgh.u(i, j) - field.fgh.u(i - 1, j)) / c_params.mesh.mesh_dx +
         (field.fgh.v(i, j) - field.fgh.v(i, j - 1)) / c_params.mesh.mesh_dy);
  }
};
} // namespace stencils
