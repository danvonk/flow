#include "sor_solver.cuh"

#include "config.hpp"
#include "fields/flow_field.hpp"
#include "fields/field.hpp"
#include "parameters.hpp"

#include <cuda.h>

#include <spdlog/spdlog.h>

struct BadInfo {
  int found;
  int i, j;
  int which;
};

__global__ void find_nan(FlowFieldView field, BadInfo *out)
{
  int nx = field.Nx, ny = field.Ny;
  int stride = field.p.Nx;

  int i = 2 + (blockIdx.x * blockDim.x + threadIdx.x);
  int j = 2 + (blockIdx.y * blockDim.y + threadIdx.y);
  if (i >= nx + 2 || j >= ny + 2)
    return;

  int c = i + stride * j;

  auto p = field.p.data[c];
  auto b = field.rhs.data[c];

  if (!isfinite(p) || !isfinite(b)) {
    if (atomicCAS(&out->found, 0, 1) == 0) {
      out->i = i;
      out->j = j;
      out->which = !isfinite(p) ? 0 : 1; // 0=p, 1=rhs
    }
  }
}

__global__ void rb_sor_pressure_step(FlowFieldView field, Real aW, Real aE,
                                     Real aS, Real aN, Real inv_aC, Real omg,
                                     int target_colour)
{
  // fluid cells only
  const auto nx = field.Nx;
  const auto ny = field.Ny;

  const auto i = 2 + (blockIdx.x * blockDim.x + threadIdx.x);
  const auto j = 2 + (blockIdx.y * blockDim.y + threadIdx.y);

  if (i >= nx + 2 || j >= ny + 2)
    return;
  if (((i + j) & 1) != target_colour)
    return;

  const auto stride = field.p.Nx;
  const auto c = i + stride * j;
  const auto w = c - 1;
  const auto e = c + 1;
  const auto s = c - stride;
  const auto n = c + stride;

  // gauss-seidel
  auto gs = inv_aC *
            (field.rhs.data[c] - aW * field.p.data[w] - aE * field.p.data[e] -
             aS * field.p.data[s] - aN * field.p.data[n]);
  auto old = field.p.data[c];
  // blended update
  field.p.data[c] = omg * gs + (1.0 - omg) * old;
}

__global__ void copy_left_right_boundary(FlowFieldView field)
{
  const auto nx = field.Nx;
  const auto ny = field.Ny;
  const auto stride = field.p.Nx;

  const auto j = 2 + (blockIdx.x * blockDim.x + threadIdx.x);
  // ignore last boundary cell
  if (j >= ny + 2)
    return;

  // two ghost cells on left
  field.p.data[stride * j] = field.p.data[2 + stride * j];
  field.p.data[stride * j + 1] = field.p.data[2 + stride * j];

  // one ghost cell on right
  field.p.data[(nx + 2) + stride * j] = field.p.data[(nx + 1) + stride * j];
}

__global__ void copy_top_bottom_boundary(FlowFieldView field)
{
  const int nx = field.Nx;
  const int ny = field.Ny;
  const int stride = field.p.Nx; // includes ghosts

  const int i = 2 + (blockIdx.x * blockDim.x + threadIdx.x);

  if (i >= nx + 2)
    return;

  // bottom side (two ghosts)
  field.p.data[i + stride * 0] = field.p.data[i + stride * 2];
  field.p.data[i + stride * 1] = field.p.data[i + stride * 2];

  // top side (one ghost)
  field.p.data[i + stride * (ny + 2)] = field.p.data[i + stride * (ny + 1)];
}

__global__ void residual_renorm(FlowFieldView field, Real aW, Real aE, Real aS,
                                Real aN, Real aC, Real *__restrict__ sumSq)
{
  const int nx = field.Nx;
  const int ny = field.Ny;
  const int stride = field.p.Nx;

  const int i = 2 + (int)(blockIdx.x * blockDim.x + threadIdx.x);
  const int j = 2 + (int)(blockIdx.y * blockDim.y + threadIdx.y);

  if (i >= nx + 2 || j >= ny + 2)
    return;

  const int c = i + stride * j;
  const int w = c - 1;
  const int e = c + 1;
  const int s = c - stride;
  const int n = c + stride;

  auto r = field.rhs.data[c] - aW * field.p.data[w] - aE * field.p.data[e] -
           aS * field.p.data[s] - aN * field.p.data[n] - aC * field.p.data[c];

  atomicAdd(sumSq, r * r);
}

int SORSolver::solve(FlowField &field)
{
                          
  const int Vx = field.cellsX();
  const int Vy = field.cellsY();
 
  const Real tol = 1e-4;
  const Real omg = 1.7;

  // field cells only
  const int nx = field.Nx();
  const int ny = field.Ny();

  const auto dx = field.params()->mesh.mesh_dx;
  const auto dy = field.params()->mesh.mesh_dy;

  const Real idx2 = 1. / (dx * dx);
  const Real idy2 = 1. / (dy * dy);

  Real aW = idx2;
  Real aE = idx2;
  Real aS = idy2;
  Real aN = idy2;
  const Real aC = -2. * (idx2 + idy2);
  const Real inv_aC = 1.0 / aC;

  dim3 block(TPB);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);

  cudaMalloc(&d_sumSq, sizeof(Real));

  dim3 gridJ((ny + TPB - 1) / TPB);
  dim3 gridI((nx + TPB - 1) / TPB);

  BadInfo *d_bad = nullptr;
  cudaMalloc(&d_bad, sizeof(BadInfo));

  int it = 0;
  for (; it < max_iters_; ++it) {
    // red
    rb_sor_pressure_step<<<grid, block>>>(field.view(), aW, aE, aS, aN, inv_aC,
                                          omg, 0);

    // black
    rb_sor_pressure_step<<<grid, block>>>(field.view(), aW, aE, aS, aN, inv_aC,
                                          omg, 1);

    copy_left_right_boundary<<<gridJ, TPB>>>(field.view());

    copy_top_bottom_boundary<<<gridI, TPB>>>(field.view());

    cudaMemset(d_sumSq, 0, sizeof(Real));
    residual_renorm<<<grid, block>>>(field.view(), aW, aE, aS, aN, aC, d_sumSq);

    Real h_sumSq = 0.0;
    cudaMemcpy(&h_sumSq, d_sumSq, sizeof(Real), cudaMemcpyDeviceToHost);

    Real resnorm = sqrt(h_sumSq / (nx * ny));
    if (resnorm <= tol) {
      cudaFree(d_sumSq);
      return it + 1;
    }
  }

  cudaFree(d_sumSq);
  return it < max_iters_ - 1;
}
