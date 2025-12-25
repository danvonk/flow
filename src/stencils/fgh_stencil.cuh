#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

#include <cmath>

#include <cuda.h>
#include <cuda_runtime.h>

namespace stencils {
struct FGHStencil {
private:
  __device__ inline Real du2dx(FlowFieldView &field, int i, int j)
  {
    const auto dx_short = 0.5 * (c_params.mesh.mesh_dx);
    const auto dx_long1 = 0.5 * (c_params.mesh.mesh_dx + c_params.mesh.mesh_dx);

    const auto u_0 = field.v.u(i, j);
    const auto u_M1 = field.v.u(i - 1, j);
    const auto u_1 = field.v.u(i + 1, 0);

    const auto kr = (u_0 + u_1) / 2.0;
    const auto kl = (u_0 + u_M1) / 2.0;

    const auto secondOrder =
        ((u_0 + u_1) * (u_0 + u_1) - (u_0 + u_M1) * (u_0 + u_M1)) /
        (4 * dx_long1);

    // Donor-cell like derivative expression. We evaluate u half-way between
    // neighboured u-components and use this as a prediction of the transport
    // direction.
    const auto firstOrder = 1.0 / (4.0 * dx_short) *
                            (kr * (u_0 + u_1) - kl * (u_M1 + u_0) +
                             fabs(kr) * (u_0 - u_1) - fabs(kl) * (u_M1 - u_0));

    // Return linear combination of central- and upwind difference
    const auto gamma = d_timestep.gamma;
    return (1.0 - gamma) * secondOrder + gamma * firstOrder;
  }

  __device__ inline Real dv2dy(FlowFieldView &field, int i, int j)
  {
    const auto dy = c_params.mesh.mesh_dy;
    const auto dyShort = 0.5 * dy;
    // const auto dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const auto dyLong1 = 0.5 * (dy + dy);

    const auto v0 = field.v.v(i, j);
    const auto vM1 = field.v.v(i, j - 1);
    const auto v1 = field.v.v(i, j + 1);

    // const auto kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort /
    // dyLong1 * v1; const auto kl = (dyLong0 - dyShort) / dyLong0 * v0 +
    // dyShort / dyLong0 * vM1;
    const auto kr = (v0 + v1) / 2;
    const auto kl = (v0 + vM1) / 2;

    /*const auto secondOrder = (((dyLong1 - dyShort) / dyLong1 * v0 +
       dyShort / dyLong1 * v1) * ((dyLong1 - dyShort) / dyLong1 * v0 + dyShort /
       dyLong1 * v1)
        - ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) *
       ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) ) / (2.0 *
       dyShort);*/

    const auto secondOrder =
        ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const auto firstOrder = 1.0 / (4.0 * dyShort) *
                            (kr * (v0 + v1) - kl * (vM1 + v0) +
                             fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const auto gamma = d_timestep.gamma;
    const auto tmp2 = (1.0 - gamma) * secondOrder + gamma * firstOrder;

    return tmp2;
  }

  __device__ inline Real duvdy(FlowFieldView &field, int i, int j)
  {

    // TODO: simplify now we only use uniform grids
    const auto hyShort =
        0.5 * c_params.mesh.mesh_dy; // Distance of corner points in
                                     // x-direction from center v-value
    const auto hyLong0 =
        0.5 * (c_params.mesh.mesh_dy +
               c_params.mesh.mesh_dy); // Distance between center
                                       // and west v-value
    const auto hyLong1 =
        0.5 * (c_params.mesh.mesh_dy +
               c_params.mesh.mesh_dy); // Distance between center and
                                       // east v-value

    // Distance of center u-value from upper edge of cell
    const auto hxShort = 0.5 * c_params.mesh.mesh_dy;
    // Distance of north and center u-value
    const auto hxLong = 0.5 * (c_params.mesh.mesh_dy + c_params.mesh.mesh_dy);

    const auto v00 = field.v.v(i, j);
    const auto v10 = field.v.v(i + 1, j);
    const auto u00 = field.v.u(i, j);
    const auto u01 = field.v.u(i, j + 1);

    const auto v0M1 = field.v.v(i, j - 1);
    const auto v1M1 = field.v.v(i + 1, j - 1);
    const auto u0M1 = field.v.u(i, j - 1);

    const auto secondOrder =
        (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10) *
             ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01) -
         ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1) *
             ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1)) /
        (2.0 * hyShort);

    const auto kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const auto kl =
        (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const auto gamma = d_timestep.gamma;
    const auto firstOrder = 1.0 / (4.0 * hyShort) *
                            (kr * (u00 + u01) - kl * (u0M1 + u00) +
                             fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const auto tmp2 = (1.0 - gamma) * secondOrder + gamma * firstOrder;

    return tmp2;
  }

  __device__ inline Real duvdx(FlowFieldView &field, int i, int j)
  {
    // TODO: simplify for uniform grids
    const auto hxShort =
        0.5 * c_params.mesh.mesh_dx; // Distance of corner points in
                                     // x-direction from center v-value
    const auto hxLong0 =
        0.5 * (c_params.mesh.mesh_dx +
               c_params.mesh.mesh_dx); // Distance between center and
                                       // west v-value
    const auto hxLong1 =
        0.5 * (c_params.mesh.mesh_dx +
               c_params.mesh.mesh_dx); // Distance between center and
                                       // east v-value

    // Distance of center u-value from upper edge of cell
    const auto hyShort = 0.5 * c_params.mesh.mesh_dy;
    // Distance of north and center u-value
    const auto hyLong = 0.5 * (c_params.mesh.mesh_dy + c_params.mesh.mesh_dy);

    const auto u00 = field.v.u(i, j);
    const auto u01 = field.v.u(i, j + 1);
    const auto v00 = field.v.v(i, j);
    const auto v10 = field.v.v(i + 1, j);

    const auto uM10 = field.v.u(i - 1, j);
    const auto uM11 = field.v.u(i - 1, j + 1);
    const auto vM10 = field.v.v(i - 1, j);

    // This a central difference expression for the first-derivative. We
    // therefore linearly interpolate u*v onto the surface of the current cell
    // (in 2D: upper left and upper right corner) and then take the central
    // difference.
    const auto secondOrder =
        (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01) *
             ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10) -
         ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11) *
             ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10)) /
        (2.0 * hxShort);

    // This is a forward-difference in donor-cell style. We apply donor cell and
    // again interpolate the velocity values (u-comp.) onto the surface of the
    // cell. We then apply the standard donor cell scheme. This will, however,
    // result in non-equal mesh spacing evaluations (in case of stretched
    // meshes).
    const auto kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
    const auto kl =
        (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

    const auto firstOrder = 1.0 / (4.0 * hxShort) *
                            (kr * (v00 + v10) - kl * (vM10 + v00) +
                             fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

    // Return linear combination of central and donor-cell difference

    const auto gamma = d_timestep.gamma;
    const auto tmp2 = (1.0 - gamma) * secondOrder + gamma * firstOrder;

    return tmp2;
  }

  __device__ inline Real d2udx2(FlowFieldView &field, int i, int j)
  {
    auto u_i = field.v.u(i, j);
    auto u_im1 = field.v.u(i - 1, j);
    auto u_ip1 = field.v.u(i + 1, j);

    auto dx = c_params.mesh.mesh_dx;

    return (u_ip1 - 2. * u_i + u_im1) / (dx * dx);
  }

  __device__ inline Real d2udy2(FlowFieldView &field, int i, int j)
  {
    auto u_i = field.v.u(i, j);
    auto u_jm1 = field.v.u(i, j - 1);
    auto u_jp1 = field.v.u(i, j + 1);

    const auto dy = c_params.mesh.mesh_dy;

    return (u_jp1 - 2. * u_i + u_jm1) / (dy * dy);
  }

  __device__ inline Real d2vdx2(FlowFieldView &field, int i, int j)
  {
    auto v_i = field.v.v(i, j);
    auto v_im1 = field.v.v(i - 1, j);
    auto v_ip1 = field.v.v(i + 1, j);

    const auto dx = c_params.mesh.mesh_dx;

    return (v_ip1 - 2. * v_i + v_im1) / (dx * dx);
  }

  __device__ inline Real d2vdy2(FlowFieldView &field, int i, int j)
  {
    auto v_i = field.v.v(i, j);
    auto v_jm1 = field.v.v(i, j - 1);
    auto v_jp1 = field.v.v(i, j + 1);

    const auto dy = c_params.mesh.mesh_dy;

    return (v_jp1 - 2. * v_i + v_jm1) / (dy * dy);
  }

public:
  __device__ void operator()(FlowFieldView &field, int i, int j)
  {
    const auto dt = d_timestep.dt;
    // compute F
    {
      auto u_current = field.v.u(i, j);
      const auto advection = -du2dx(field, i, j) - duvdy(field, i, j);
      const auto force = c_params.sim.gx;
      auto diffusion = d2udx2(field, i, j) + d2udy2(field, i, j);
      diffusion *= (1.0 / c_params.sim.reynolds);

      field.fgh.u(i, j) = u_current + dt * (advection + diffusion + force);
    }

    // compute G
    {
      auto v_current = field.v.v(i, j);
      const auto advection = -dv2dy(field, i, j) - duvdx(field, i, j);
      const auto force = c_params.sim.gy;
      auto diffusion = d2vdx2(field, i, j) + d2vdy2(field, i, j);
      diffusion *= (1.0 / c_params.sim.reynolds);

      field.fgh.v(i, j) = v_current + dt * (advection + diffusion + force);
    }
  }
};

} // namespace stencils
