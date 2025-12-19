#pragma once

#include <memory>

#include "fields/flow_field.hpp"

class SimulationImpl;

class Simulation {
public:
  Simulation(FlowField &field);
  ~Simulation();

  void run();

private:
  SimulationImpl *impl_;
};
