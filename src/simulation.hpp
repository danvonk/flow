#pragma once

class Simulation {
public:
  Simulation();

  void initFlowField();
  void solveTimestep();
  void setTimestep();
};
