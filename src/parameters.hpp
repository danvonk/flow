#pragma once

#include "config.hpp"

namespace params {
struct Timestep {
  Real dt = 0;
  Real tau = 0;
};

struct Simulation {
  Real final_time = 0;
  Real reynolds = 0;
};

} // namespace params

class Parameters {};
