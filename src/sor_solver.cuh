#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

class SORSolver {
public:
  int solve(FlowField &field);

  inline auto max_iters() const { return max_iters_; };

private:
  Real *d_sumSq = nullptr;
  const int max_iters_ = 100000;
};
