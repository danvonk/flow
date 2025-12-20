#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

class SORSolver {
public:
  SORSolver(FlowField &field);
  void solve();

private:
  FlowField &field_;
};
