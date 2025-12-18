#pragma once

#include <memory>

#include "fields/flow_field.hpp"

class SimulationImpl;

class Simulation {
public:
  Simulation(FlowField &field);
  ~Simulation();

  void initFlowField();
  void solveTimestep();
  void setTimestep();

private:
  SimulationImpl *impl_;
};
