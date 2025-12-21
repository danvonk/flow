#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"

class SORSolver {
public:
  int solve(FlowFieldView field);
};
