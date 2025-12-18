#pragma once

#include "config.hpp"

#include "fields/flow_field.hpp"
#include "stencil.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

namespace iter_kernels {
template <typename Stencil>
__global__ void apply_stencil_kernel(FlowFieldView flow, Stencil stencil,
                                     int i0, int i1, int j0, int j1)
{
  int i = i0 + blockIdx.x * blockDim.x + threadIdx.x;
  int j = j0 + blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= i1 || j >= j1)
    return;

  stencil(flow, i, j); // call stencil on device
}

template <typename Stencil>
__global__ void apply_left_wall(FlowFieldView flow, Stencil stencil, int i_wall,
                                int j0, int j1)
{
  int j = j0 + blockIdx.x * blockDim.x + threadIdx.x;
  if (j < j1)
    stencil.apply_left_wall(flow, i_wall, j);
}

template <typename Stencil>
__global__ void apply_right_wall(FlowFieldView flow, Stencil stencil,
                                 int i_wall, int j0, int j1)
{
  int j = j0 + blockIdx.x * blockDim.x + threadIdx.x;
  if (j < j1)
    stencil.apply_right_wall(flow, i_wall, j);
}

template <typename Stencil>
__global__ void apply_bottom_wall(FlowFieldView flow, Stencil stencil, int i0,
                                  int i1, int j_wall)
{
  int i = i0 + blockIdx.x * blockDim.x + threadIdx.x;
  if (i < i1)
    stencil.apply_bottom_wall(flow, i, j_wall);
}

template <typename Stencil>
__global__ void apply_top_wall(FlowFieldView flow, Stencil stencil, int i0,
                               int i1, int j_wall)
{
  int i = i0 + blockIdx.x * blockDim.x + threadIdx.x;
  if (i < i1)
    stencil.apply_top_wall(flow, i, j_wall);
}
} // namespace iter_kernels

template <FieldStencil2D Stencil> class FieldIterator {
public:
  FieldIterator(FlowFieldView field, Stencil stencil, int low_offset = 0,
                int high_offset = 0)
      : field_(field),
        stencil_(stencil),
        low_offset_(low_offset),
        high_offset_(high_offset)
  {
  }

  void iterate(cudaStream_t stream = 0)
  {
    int cellsX = field_.cellsX;
    int cellsY = field_.cellsY;

    const int i0 = 1 + low_offset_;
    const int i1 = (cellsX - 1) + high_offset_;
    const int j0 = 1 + low_offset_;
    const int j1 = (cellsY - 1) + high_offset_;

    dim3 block(16, 16, 1);
    dim3 grid(((i1 - i0) + block.x - 1) / block.x,
              ((j1 - j0) + block.y - 1) / block.y, 1);

    iter_kernels::apply_stencil_kernel<<<grid, block, 0, stream>>>(
        field_, stencil_, i0, i1, j0, j1);
  }

private:
  FlowFieldView field_;
  Stencil stencil_;
  int low_offset_;
  int high_offset_;
};

template <BoundaryStencil2D Stencil> class BoundaryIterator {
public:
  BoundaryIterator(FlowFieldView field, Stencil stencil, int low_offset = 0,
                   int high_offset = 0)
      : field_(field),
        stencil_(stencil),
        low_offset_(low_offset),
        high_offset_(high_offset)
  {
  }

  void iterate(cudaStream_t stream = 0)
  {
    int cellsX = field_.cellsX;
    int cellsY = field_.cellsY;

    int j0 = low_offset_;
    int j1 = cellsY + high_offset_;
    int i0 = low_offset_;
    int i1 = cellsX + high_offset_;

    int left_i = low_offset_;
    int right_i = cellsX + high_offset_ - 1;
    int bot_j = low_offset_;
    int top_j = cellsY + high_offset_ - 1;

    constexpr int TPB = 256; // threads per block

    int nJ = j1 - j0;
    dim3 blockJ(TPB);
    dim3 gridJ((nJ + TPB - 1) / TPB);

    iter_kernels::apply_left_wall<<<gridJ, blockJ, 0, stream>>>(
        field_, stencil_, left_i, j0, j1);
    iter_kernels::apply_right_wall<<<gridJ, blockJ, 0, stream>>>(
        field_, stencil_, right_i, j0, j1);

    int nI = i1 - i0;
    dim3 blockI(TPB);
    dim3 gridI((nI + TPB - 1) / TPB);

    iter_kernels::apply_bottom_wall<<<gridI, blockI, 0, stream>>>(
        field_, stencil_, i0, i1, bot_j);
    iter_kernels::apply_top_wall<<<gridI, blockI, 0, stream>>>(field_, stencil_,
                                                               i0, i1, top_j);
  }

private:
  FlowFieldView field_;
  Stencil stencil_;
  int low_offset_;
  int high_offset_;
};
