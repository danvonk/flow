#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

class SORSolver {
public:
  int solve(FlowField &field);

private:
  Real *d_sumSq = nullptr;
};
