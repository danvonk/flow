#pragma once

#include "config.hpp"

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

struct Timestep {
  Real dt = 1.0;
  Real tau = 1.0;
  Real gamma = 1.0; // central vs upwind scheme balance coeff
};

struct Simulation {
  Real final_time = 100;
  Real reynolds = 1000;

  // gravity
  Real gx = 0.0;
  Real gy = 0.0;
};

struct Mesh {
  Real sizeX = 1.;
  Real sizeY = 1.;

  int cellsX = 100;
  int cellsY = 100;

  Real mesh_dx;
  Real mesh_dy;
};

struct Wall {
  Real scalar_left = 0.0;
  Real scalar_right = 0.0;
  Real scalar_top = 0.0;
  Real scalar_bottom = 0.0;

  // define the velocities at the boundaries
  Real vector_left[2] = {0, 0};
  Real vector_right[2] = {0, 0};
  Real vector_top[2] = {1, 0};
  Real vector_bottom[2] = {0, 0};

  BoundaryType boundary_left = BoundaryType::DIRICHLET;
  BoundaryType boundary_right = BoundaryType::DIRICHLET;
  BoundaryType boundary_top = BoundaryType::DIRICHLET;
  BoundaryType boundary_bottom = BoundaryType::DIRICHLET;
};

struct Parameters {
  Timestep timestep;
  Simulation sim;
  Mesh mesh;
  Wall wall;
};

Parameters load_config(const std::string &file_path);

} // namespace params
