#include "cfl_helpers.cuh"
#include "simulation_impl.cuh"
#include "vtk_writer.hpp"

#include <spdlog/spdlog.h>

SimulationImpl::SimulationImpl(FlowField &field)
    : flow_(field),
      fgh_iterator_(field.view(), stencils::FGHStencil()),
      wall_v_iterator_(field.view(), stencils::MovingWallVelocityStencil()),
      wall_fgh_iterator_(field.view(), stencils::MovingWallFGHStencil()),
      rhs_iterator_(field.view(), stencils::RHSStencil()),
      velocity_stencil_(field.view(), stencils::VelocityStencil()),
      obstacle_stencil_(field.view(), stencils::ObstacleStencil()),
      solver_(SORSolver(field))
{
}

void SimulationImpl::solveTimestep()
{
  fgh_iterator_.iterate();
  wall_fgh_iterator_.iterate();
  rhs_iterator_.iterate();

  // TODO: solve for pressure
  solver_.solve();

  // compute velocity
  velocity_stencil_.iterate();
  obstacle_stencil_.iterate();

  // iterate for velocities on boundary
  wall_v_iterator_.iterate();
}

void SimulationImpl::run()
{
  Real time = 0.0;
  int timesteps = 0;
  VTKWriter writer;

  while (time < flow_.params()->sim.final_time) {
    auto new_dt = calc_new_timestep(flow_.view());
    solveTimestep();

    time += new_dt;
    timesteps++;

    writer.write_flow_field(flow_, time);
    spdlog::info("Current time: {}\t Timestep: {}", time, new_dt);
  }
}
