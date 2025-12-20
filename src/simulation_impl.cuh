#pragma once

#include "fields/flow_field.hpp"
#include "parameters.hpp"
#include "sor_solver.cuh"
#include "stencils/fgh_stencil.cuh"
#include "stencils/iterator.cuh"
#include "stencils/moving_wall_stencil.cuh"
#include "stencils/obstacle_stencil.cuh"
#include "stencils/rhs_stencil.cuh"
#include "stencils/velocity_stencil.cuh"

class SimulationImpl {
public:
  SimulationImpl(FlowField &field);

  void run();

private:
  void initFlowField() {};
  void solveTimestep();

private:
  FlowField &flow_;

  FieldIterator<stencils::FGHStencil> fgh_iterator_;
  BoundaryIterator<stencils::MovingWallVelocityStencil> wall_v_iterator_;
  BoundaryIterator<stencils::MovingWallFGHStencil> wall_fgh_iterator_;
  FieldIterator<stencils::RHSStencil> rhs_iterator_;
  FieldIterator<stencils::VelocityStencil> velocity_stencil_;
  FieldIterator<stencils::ObstacleStencil> obstacle_stencil_;

  SORSolver solver_;
};
