#include "simulation_impl.cuh"

SimulationImpl::SimulationImpl(FlowField &field)
    : flow_(field),
      fgh_iterator_(field.view(), FGHStencil()),
      wall_v_iterator_(field.view(), MovingWallVelocityStencil()),
      wall_fgh_iterator_(field.view(), MovingWallFGHStencil())
{
}

void SimulationImpl::setTimestep() {}

void SimulationImpl::solveTimestep()
{
  setTimestep();
  fgh_iterator_.iterate();
  wall_fgh_iterator_.iterate();
}
