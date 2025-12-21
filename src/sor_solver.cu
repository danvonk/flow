#include "sor_solver.cuh"

#include "parameters.hpp"

#include <cuda.h>

__global__ void rb_sor_pressure_step(FlowFieldView field, Real aW, Real aE,
                                     Real aS, Real aN, Real inv_aC, Real omg,
                                     int target_colour)
{
  const auto nx = field.Nx;
  const auto ny = field.Ny;

  const auto i = 2 + (blockIdx.x * blockDim.x + threadIdx.x);
  const auto j = 2 + (blockIdx.y * blockDim.y + threadIdx.y);

  if (i > nx + 1 || j > ny + 1)
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
  // update
  field.p.data[c] = omg * gs * (1. - omg) * old;
}

__global__ void copy_left_right_boundary(FlowFieldView field)
{
  const auto nx = field.Nx;
  const auto ny = field.Ny;
  const auto stride = field.p.Nx;

  const auto j = 2 + (blockIdx.x + blockDim.x + threadIdx.x);
  // ignore last boundary cell
  if (j > ny + 1)
    return;

  // two ghost cells on left
  field.p.data[stride * j] = field.p.data[2 + stride * j];
  field.p.data[stride * j + 1] = field.p.data[2 + stride * j];

  // one ghost cell on right
  field.p.data[(nx + 2) + stride * j] = field.p.data[(nx + 1) + stride * j];
}

__global__ void copy_top_bottom_boundary(FlowFieldView field) {}

__global__ void residual_renorm(FlowFieldView field, Real aW, Real aE, Real aS,
                                Real aN, Real aC)
{
}

// int SORSolver::solve(FlowFieldView field) {}

int SORSolver::solve(FlowFieldView field)
{

  Real resnorm = std::numeric_limits<Real>::max();
  Real tol = 1e-4;

  double omg = 1.7;
  int iterations = -1;
  int it = 0;

  // only field cells
  int nx = field.Nx;
  int ny = field.Ny;

  const auto dx = c_params.mesh.mesh_dx;
  const auto dy = c_params.mesh.mesh_dy;

  do {
    for (int j = 2; j < ny + 2; j++) {
      for (int i = 2; i < nx + 2; i++) {

        const Real dx_W = 0.5 * (dx + dx);
        const Real dx_E = 0.5 * (dx + dx);
        const Real dx_S = 0.5 * (dy + dy);
        const Real dx_N = 0.5 * (dy + dy);

        const Real a_W = 2.0 / (dx_W * (dx_W + dx_E));
        const Real a_E = 2.0 / (dx_E * (dx_W + dx_E));
        const Real a_N = 2.0 / (dx_N * (dx_N + dx_S));
        const Real a_S = 2.0 / (dx_S * (dx_N + dx_S));
        const Real a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);

        const Real gaussSeidel =
            1.0 / a_C *
            (field.rhs.get(i, j) - a_W * field.p.get(i - 1, j) -
             a_E * field.p.get(i + 1, j) - a_S * field.p.get(i, j - 1) -
             a_N * field.p.get(i, j + 1));

        field.p.get(i, j) = omg * gaussSeidel + (1.0 - omg) * field.p.get(i, j);
      }
    }

    resnorm = 0.0;
    for (int j = 2; j < ny + 2; j++) {
      for (int i = 2; i < nx + 2; i++) {

        const Real dx_W = 0.5 * (dx + dx);
        const Real dx_E = 0.5 * (dx + dx);
        const Real dx_S = 0.5 * (dy + dy);
        const Real dx_N = 0.5 * (dy + dy);

        const Real a_W = 2.0 / (dx_W * (dx_W + dx_E));
        const Real a_E = 2.0 / (dx_E * (dx_W + dx_E));
        const Real a_N = 2.0 / (dx_N * (dx_N + dx_S));
        const Real a_S = 2.0 / (dx_S * (dx_N + dx_S));
        const Real a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);

        const Real residual =
            field.rhs.get(i, j) - a_W * field.p.get(i - 1, j) -
            a_E * field.p.get(i + 1, j) - a_S * field.p.get(i, j - 1) -
            a_N * field.p.get(i, j + 1) - a_C * field.p.get(i, j);
        resnorm += residual * residual;
      }
    }
    resnorm = sqrt(resnorm / (nx * ny));

    for (int j = 2; j < ny + 2; j++) {
      field.p.get(1, j) = field.p.get(2, j);
      field.p.get(nx + 2, j) = field.p.get(nx + 1, j);
    }

    for (int i = 2; i < nx + 2; i++) {
      field.p.get(i, 1) = field.p.get(i, 2);
      field.p.get(i, ny + 2) = field.p.get(i, ny + 1);
    }
    iterations--;
    it++;
  } while (resnorm > tol && iterations);

  return it;
}
