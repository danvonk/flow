#include "simulation_impl.cuh"
#include "vtk_writer.hpp"

#include <spdlog/spdlog.h>

SimulationImpl::SimulationImpl(FlowField &field)
    : flow_(field),
      fgh_iterator_(field.view(), stencils::FGHStencil()),
      wall_v_iterator_(field.view(), stencils::MovingWallVelocityStencil()),
      wall_fgh_iterator_(field.view(), stencils::MovingWallFGHStencil()),
      rhs_iterator_(field.view(), stencils::RHSStencil())
{
}

void SimulationImpl::setTimestep() {}

void SimulationImpl::solveTimestep()
{
  spdlog::info("Solving timestep...");
  setTimestep();
  fgh_iterator_.iterate();
  wall_fgh_iterator_.iterate();
  rhs_iterator_.iterate();

  // TODO: solve for pressure
}

void SimulationImpl::run()
{
  Real time = 0.0;
  solveTimestep();

  VTKWriter writer;
  writer.write_flow_field(flow_, 0.0);
}
