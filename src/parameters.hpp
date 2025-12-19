#pragma once

#include "config.hpp"

#include <cuda.h>
#include <string>

namespace params {
// Types of boundary
enum class BoundaryType {
  DIRICHLET,
  // Only need this for other scenarios...
  PERIODIC,
  PARALLEL_BOUNDARY,
  NEUMANN
};

static constexpr int OBSTACLE_SELF = 1 << 0;
static constexpr int OBSTACLE_LEFT = 1 << 1;
static constexpr int OBSTACLE_RIGHT = 1 << 2;
static constexpr int OBSTACLE_BOTTOM = 1 << 3;
static constexpr int OBSTACLE_TOP = 1 << 4;
static constexpr int OBSTACLE_FRONT = 1 << 5;
static constexpr int OBSTACLE_BACK = 1 << 6;

struct Simulation {
  Real final_time;
  Real reynolds;

  // gravity
  Real gx;
  Real gy;
};

struct Mesh {
  // physical dimensions
  Real lengthX;
  Real lengthY;

  // number of (non-ghost) cells
  int cellsX;
  int cellsY;

  Real mesh_dx;
  Real mesh_dy;
};

struct Wall {
  Real scalar_left;
  Real scalar_right;
  Real scalar_top;
  Real scalar_bottom;

  // define the velocities at the boundaries
  Real vector_left[2];
  Real vector_right[2];
  Real vector_top[2];
  Real vector_bottom[2];

  BoundaryType boundary_left;
  BoundaryType boundary_right;
  BoundaryType boundary_top;
  BoundaryType boundary_bottom;
};

struct Parameters {
  Simulation sim;
  Mesh mesh;
  Wall wall;
};

struct Timestep {
  Real dt;
  Real tau;
  Real gamma; // central vs upwind scheme balance coeff
};

Parameters load_config(const std::string &file_path);

} // namespace params

// These live on the CUDA device
extern __constant__ params::Parameters c_params;
extern __device__ params::Timestep d_timestep;
