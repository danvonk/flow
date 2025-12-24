#include "sor_solver.cuh"

#include "config.hpp"
#include "fields/flow_field.hpp"
#include "fields/field.hpp"
#include "parameters.hpp"

#include <cuda.h>
#include <vector>
#include <iostream>

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

  // iterate interior i indices: i=2..nx+1 (count nx)
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

  atomicAdd(sumSq, (Real)r * (Real)r);
}

int SORSolver::solve(FlowField &field)
{

  std::vector<Real> rhs_data;
  std::vector<Real> velocity_data;
  field.rhs()->to_host(rhs_data);
  field.fgh()->to_host(velocity_data);

  auto nxtot = field.rhs()->Nx();
  auto nytot = field.rhs()->Ny();
  // spdlog::info("RHS field {}x{}", nxtot, nytot);
  
  // for (int j = 0; j < nytot; j++) {
  //   for (int i = 0; i < nxtot; i++) {
  //       std::cout << "(" << i << "," << j << "): " << rhs_data[i + nxtot * j] << " ";
  //   }
  //   std::cout << '\n';
  // }
  // spdlog::info("===");
                           
  const int Vx = field.cellsX();
  const int Vy = field.cellsY();

  // spdlog::info("FGH field {}x{}", Vx, Vy);
  // for (int j = 0; j < Vy; ++j) {
  //   for (int i = 0; i < Vx; ++i) {
  //     auto here = 2 * (i + Vx * j);
  //     std::cout << velocity_data[here + 0] << "," << velocity_data[here + 1]
  //               << " ";
  //   }
  //   std::cout << '\n';
  // }

  
  const int maxIters = 1000;
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

  dim3 block(16, 16);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);

  cudaMalloc(&d_sumSq, sizeof(Real));

  dim3 gridJ((ny + TPB - 1) / TPB);
  dim3 gridI((nx + TPB - 1) / TPB);

  BadInfo *d_bad = nullptr;
  cudaMalloc(&d_bad, sizeof(BadInfo));

  for (int it = 0; it < maxIters; ++it) {

    find_nan<<<grid, block>>>(field.view(), d_bad);
    BadInfo h_bad{};
    cudaMemcpy(&h_bad, d_bad, sizeof(BadInfo), cudaMemcpyDeviceToHost);
    if (h_bad.found) {
      const char *what = (h_bad.which == 0) ? "pressure" : "rhs";
      spdlog::error("NaN found in {} at i={}, j={}", what, h_bad.i, h_bad.j);
    }
    cudaMemset(d_bad, 0, sizeof(BadInfo));

    // red
    rb_sor_pressure_step<<<grid, block>>>(field.view(), aW, aE, aS, aN, inv_aC,
                                          omg, 0);
    cudaDeviceSynchronize();
    auto err = cudaGetLastError();
    if (err != cudaSuccess)
      printf("CUDA error: %s\n", cudaGetErrorString(err));

    // black
    rb_sor_pressure_step<<<grid, block>>>(field.view(), aW, aE, aS, aN, inv_aC,
                                          omg, 1);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess)
      printf("CUDA error: %s\n", cudaGetErrorString(err));

    copy_left_right_boundary<<<gridJ, TPB>>>(field.view());
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess)
      printf("CUDA error: %s\n", cudaGetErrorString(err));

    copy_top_bottom_boundary<<<gridI, TPB>>>(field.view());
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess)
      printf("CUDA error: %s\n", cudaGetErrorString(err));

    cudaMemset(d_sumSq, 0, sizeof(Real));
    residual_renorm<<<grid, block>>>(field.view(), aW, aE, aS, aN, aC, d_sumSq);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess)
      printf("CUDA error: %s\n", cudaGetErrorString(err));

    Real h_sumSq = 0.0;
    cudaMemcpy(&h_sumSq, d_sumSq, sizeof(Real), cudaMemcpyDeviceToHost);

    Real resnorm = sqrt(h_sumSq / (nx * ny));
    if (resnorm <= tol) {
      cudaFree(d_sumSq);
      return it + 1;
    }
  }

  cudaFree(d_sumSq);
  return maxIters;
}

// int SORSolver::solve(FlowFieldView field)
// {

//   Real resnorm = std::numeric_limits<Real>::max();
//   Real tol = 1e-4;

//   Real omg = 1.7;
//   int iterations = -1;
//   int it = 0;

//   // only field cells
//   int nx = field.Nx;
//   int ny = field.Ny;

//   const auto dx = c_params.mesh.mesh_dx;
//   const auto dy = c_params.mesh.mesh_dy;

//   do {
//     for (int j = 2; j < ny + 2; j++) {
//       for (int i = 2; i < nx + 2; i++) {

//         const Real dx_W = 0.5 * (dx + dx);
//         const Real dx_E = 0.5 * (dx + dx);
//         const Real dx_S = 0.5 * (dy + dy);
//         const Real dx_N = 0.5 * (dy + dy);

//         const Real a_W = 2.0 / (dx_W * (dx_W + dx_E));
//         const Real a_E = 2.0 / (dx_E * (dx_W + dx_E));
//         const Real a_N = 2.0 / (dx_N * (dx_N + dx_S));
//         const Real a_S = 2.0 / (dx_S * (dx_N + dx_S));
//         const Real a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);

//         const Real gaussSeidel =
//             1.0 / a_C *
//             (field.rhs.get(i, j) - a_W * field.p.get(i - 1, j) -
//              a_E * field.p.get(i + 1, j) - a_S * field.p.get(i, j - 1) -
//              a_N * field.p.get(i, j + 1));

//         field.p.get(i, j) = omg * gaussSeidel + (1.0 - omg) * field.p.get(i,
//         j);
//       }
//     }

//     resnorm = 0.0;
//     for (int j = 2; j < ny + 2; j++) {
//       for (int i = 2; i < nx + 2; i++) {

//         const Real dx_W = 0.5 * (dx + dx);
//         const Real dx_E = 0.5 * (dx + dx);
//         const Real dx_S = 0.5 * (dy + dy);
//         const Real dx_N = 0.5 * (dy + dy);

//         const Real a_W = 2.0 / (dx_W * (dx_W + dx_E));
//         const Real a_E = 2.0 / (dx_E * (dx_W + dx_E));
//         const Real a_N = 2.0 / (dx_N * (dx_N + dx_S));
//         const Real a_S = 2.0 / (dx_S * (dx_N + dx_S));
//         const Real a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);

//         const Real residual =
//             field.rhs.get(i, j) - a_W * field.p.get(i - 1, j) -
//             a_E * field.p.get(i + 1, j) - a_S * field.p.get(i, j - 1) -
//             a_N * field.p.get(i, j + 1) - a_C * field.p.get(i, j);
//         resnorm += residual * residual;
//       }
//     }
//     resnorm = sqrt(resnorm / (nx * ny));

//     for (int j = 2; j < ny + 2; j++) {
//       field.p.get(1, j) = field.p.get(2, j);
//       field.p.get(nx + 2, j) = field.p.get(nx + 1, j);
//     }

//     for (int i = 2; i < nx + 2; i++) {
//       field.p.get(i, 1) = field.p.get(i, 2);
//       field.p.get(i, ny + 2) = field.p.get(i, ny + 1);
//     }
//     iterations--;
//     it++;
//   } while (resnorm > tol && iterations);

//   return it;
// }
