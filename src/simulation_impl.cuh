#pragma once

#include "fields/flow_field.hpp"

#include "stencils/fgh_stencil.cuh"
#include "stencils/iterator.cuh"
#include "stencils/moving_wall_stencil.cuh"

class SimulationImpl {
public:
  SimulationImpl(FlowField &field);
  void initFlowField() {};

  void solveTimestep();
  void setTimestep();

protected:
  FlowField &flow_;
  FieldIterator<FGHStencil> fgh_iterator_;
  BoundaryIterator<MovingWallVelocityStencil> wall_v_iterator_;
  BoundaryIterator<MovingWallFGHStencil> wall_fgh_iterator_;
};
