#pragma once

#include "fields/flow_field.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

namespace stencils {
struct MovingWallVelocityStencil {

  __device__ void apply_left_wall(FlowFieldView &field, int i, int j)
  {
    field.v.u(i, j) = field.params.wall.vector_left[0];
    field.v.v(i, j) =
        2 * field.params.wall.vector_left[1] - field.v.v(i + 1, j);
  }

  __device__ void apply_right_wall(FlowFieldView &field, int i, int j)
  {
    field.v.u(i - 1, j) = field.params.wall.vector_right[0];
    field.v.v(i, j) =
        2 * field.params.wall.vector_right[1] - field.v.v(i - 1, j);
  }

  __device__ void apply_top_wall(FlowFieldView &field, int i, int j)
  {
    field.v.u(i, j) = 2 * field.params.wall.vector_top[0] - field.v.u(i, j - 1);
    field.v.v(i, j - 1) = field.params.wall.vector_top[1];
  }

  __device__ void apply_bottom_wall(FlowFieldView &field, int i, int j)
  {
    field.v.u(i, j) =
        2 * field.params.wall.vector_bottom[0] - field.v.u(i, j + 1);
    field.v.v(i, j) = field.params.wall.vector_bottom[1];
  }
};

struct MovingWallFGHStencil {

  __device__ void apply_left_wall(FlowFieldView &field, int i, int j)
  {
    field.fgh.u(i, j) = field.params.wall.vector_left[0];
  }

  __device__ void apply_right_wall(FlowFieldView &field, int i, int j)
  {
    field.fgh.u(i - 1, j) = field.params.wall.vector_right[0];
  }

  __device__ void apply_bottom_wall(FlowFieldView &field, int i, int j)
  {
    field.fgh.v(i, j) = field.params.wall.vector_bottom[1];
  }

  __device__ void apply_top_wall(FlowFieldView &field, int i, int j)
  {
    field.fgh.v(i, j - 1) = field.params.wall.vector_top[1];
  }
};
} // namespace stencils
