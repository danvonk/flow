#include "cfl_helpers.cuh"
#include "simulation_impl.cuh"
#include "vtk_writer.hpp"

#include <chrono>
#include <spdlog/spdlog.h>

SimulationImpl::SimulationImpl(FlowField &field)
    : flow_(field),
      fgh_iterator_(field.view(), stencils::FGHStencil()),
      wall_v_iterator_(field.view(), stencils::MovingWallVelocityStencil()),
      wall_fgh_iterator_(field.view(), stencils::MovingWallFGHStencil()),
      rhs_iterator_(field.view(), stencils::RHSStencil()),
      velocity_stencil_(field.view(), stencils::VelocityStencil()),
      obstacle_stencil_(field.view(), stencils::ObstacleStencil())
{
}

void SimulationImpl::solveTimestep()
{
  fgh_iterator_.iterate();
  wall_fgh_iterator_.iterate();
  rhs_iterator_.iterate();

  // TODO: solve for pressure
  const int iters = solver_.solve(flow_.view());
  spdlog::debug("SORSolver needed {} iterations", iters);

  // compute velocity
  velocity_stencil_.iterate();
  obstacle_stencil_.iterate();

  // iterate for velocities on boundary
  wall_v_iterator_.iterate();
}

void SimulationImpl::run()
{
  Real time = 0.0;
  Real time_vtk = flow_.params()->sim.vtk_interval;
  int timesteps = 0;
  VTKWriter writer;

  auto start = std::chrono::steady_clock::now();
  while (time < flow_.params()->sim.final_time) {
    auto new_dt = calc_new_timestep(flow_.view());
    solveTimestep();

    time += new_dt;
    timesteps++;

    if (time_vtk <= time) {
      writer.write_flow_field(flow_, timesteps);
      time_vtk += flow_.params()->sim.vtk_interval;

      spdlog::info("Current time: {}\t Timestep: {}", time, new_dt);
    }
  }

  auto end = std::chrono::steady_clock::now();
  spdlog::info(
      "Finished simulation with a duration of {}ms",
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count());
}
