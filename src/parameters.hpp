#pragma once

#include "config.hpp"

#include <string>

namespace params {
struct Timestep {
  Real dt = 1.0;
  Real tau = 1.0;
};

struct Simulation {
  Real final_time = 100;
  Real reynolds = 1000;
};

struct Mesh {
  Real sizeX = 1;
  Real sizeY = 1;
  Real cellsX = 100;
  Real cellsY = 100;
};

struct Parameters {
  Timestep timestep;
  Simulation sim;
  Mesh mesh;

  void load_config(const std::string &file_path);
};

} // namespace params
