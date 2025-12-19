#include "config.hpp"
#include "parameters.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

using params::Parameters;

__device__ params::Timestep d_timestep;
__constant__ params::Parameters c_params;
static_assert(std::is_trivially_copyable_v<Parameters>);
static_assert(std::is_standard_layout_v<Parameters>);

Parameters params::load_config(const std::string &file_path)
{
  if (file_path.empty())
    return {};

  pt::ptree tree;
  pt::read_ini(file_path, tree);

  Parameters p;
  Timestep ts;

  ts.dt = tree.get<Real>("timestep.dt", 1.0);
  ts.tau = tree.get<Real>("timestep.tau", 1.0);
  ts.gamma = tree.get<Real>("timestep.gamma", 1.0);

  // Sim params
  p.sim.final_time = tree.get<Real>("simulation.final_time", 100);
  p.sim.reynolds = tree.get<Real>("simulation.reynolds", 1000);
  p.sim.gx = tree.get<Real>("simulation.gx", 0.0);
  p.sim.gy = tree.get<Real>("simulation.gy", 0.0);

  // Mesh params
  p.mesh.sizeX = tree.get<Real>("mesh.lengthX", 1);
  p.mesh.sizeY = tree.get<Real>("mesh.lengthY", 1);
  p.mesh.cellsX = tree.get<int>("mesh.cellsX", 100);
  p.mesh.cellsY = tree.get<int>("mesh.cellsY", 100);
  p.mesh.mesh_dx = p.mesh.sizeX / p.mesh.cellsX;
  p.mesh.mesh_dy = p.mesh.sizeY / p.mesh.cellsY;

  // Wall params
  p.wall.scalar_left = tree.get<Real>("wall.scalar_left", 0.0);
  p.wall.scalar_right = tree.get<Real>("wall.scalar_right", 0.0);
  p.wall.scalar_top = tree.get<Real>("wall.scalar_top", 0.0);
  p.wall.scalar_bottom = tree.get<Real>("wall.scalar_bottom", 0.0);

  p.wall.vector_left[0] = 0.0;
  p.wall.vector_left[1] = 0.0;
  p.wall.vector_right[0] = 0.0;
  p.wall.vector_right[1] = 0.0;
  p.wall.vector_top[0] = 1.0;
  p.wall.vector_top[1] = 0.0;
  p.wall.vector_bottom[0] = 0.0;
  p.wall.vector_bottom[1] = 0.0;

  p.wall.boundary_left = BoundaryType::DIRICHLET;
  p.wall.boundary_right = BoundaryType::DIRICHLET;
  p.wall.boundary_top = BoundaryType::DIRICHLET;
  p.wall.boundary_bottom = BoundaryType::DIRICHLET;

  // copy to device symbol
  cudaMemcpyToSymbol(d_timestep, &ts, sizeof(ts));
  cudaMemcpyToSymbol(c_params, &p, sizeof(p));

  return p;
}
