#include "cfl_helpers.cuh"

#include "config.hpp"
#include "parameters.hpp"

#include <driver_types.h>
#include <limits>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>

#include <cmath>
#include <cuda_runtime.h>

struct MaxUFunctor {
  const Real *data;
  __host__ __device__ Real operator()(int idx) const
  {
    return fabs(data[2 * idx + 0]);
  }
};

struct MaxVFunctor {
  const Real *data;
  __host__ __device__ Real operator()(int idx) const
  {
    return fabs(data[2 * idx + 1]);
  }
};

__global__ void update_dt_kernel(double u_max, double v_max, double eps)
{

  const auto dt_u = c_params.mesh.mesh_dx / (u_max + eps);
  const auto dt_v = c_params.mesh.mesh_dy / (v_max + eps);

  const auto factor = 1. / (c_params.mesh.mesh_dx * c_params.mesh.mesh_dx) +
                      1. / (c_params.mesh.mesh_dy * c_params.mesh.mesh_dy);
  const auto dt_nu = c_params.sim.reynolds / (2 * factor);

  auto dt = fmin(fmin(dt_u, dt_v), dt_nu);
  dt *= d_timestep.tau;

  d_timestep.dt = dt;
}

Real calc_new_timestep(FlowFieldView field, cudaStream_t stream)
{
  const auto N = field.v.Nx * field.v.Ny;

  auto begin = thrust::make_counting_iterator<int>(0);
  auto end = begin + N;

  auto u_max = thrust::transform_reduce(begin, end, MaxUFunctor{field.v.data},
                                        0., thrust::maximum<Real>());
  auto v_max = thrust::transform_reduce(begin, end, MaxVFunctor{field.v.data},
                                        0., thrust::maximum<Real>());

  update_dt_kernel<<<1, 1, 0, stream>>>(u_max, v_max,
                                        std::numeric_limits<Real>::epsilon());

  Real new_dt = 0.;
  cudaMemcpyFromSymbol(&new_dt, d_timestep.dt, sizeof(Real));

  return new_dt;
}
