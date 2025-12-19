#pragma once

#include "fields/flow_field.hpp"

#include "stencils/fgh_stencil.cuh"
#include "stencils/iterator.cuh"
#include "stencils/moving_wall_stencil.cuh"
#include "stencils/rhs_stencil.cuh"

#include "vtk_writer.hpp"

class SimulationImpl {
public:
  SimulationImpl(FlowField &field);

  void run();

private:
  void initFlowField() {};
  void solveTimestep();
  void setTimestep();

protected:
  FlowField &flow_;
  FieldIterator<stencils::FGHStencil> fgh_iterator_;
  BoundaryIterator<stencils::MovingWallVelocityStencil> wall_v_iterator_;
  BoundaryIterator<stencils::MovingWallFGHStencil> wall_fgh_iterator_;
  FieldIterator<stencils::RHSStencil> rhs_iterator_;
};
