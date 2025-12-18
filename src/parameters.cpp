#include "parameters.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

using params::Parameters;

Parameters params::load_config(const std::string &file_path)
{
  if (file_path.empty())
    return {};

  pt::ptree tree;
  pt::read_ini(file_path, tree);

  Parameters p;

  p.timestep.dt = tree.get<Real>("timestep.dt", 1.0);
  p.timestep.tau = tree.get<Real>("timestep.tau", 1.0);

  p.mesh.sizeX = tree.get<Real>("mesh.lengthX", 1);
  p.mesh.sizeY = tree.get<Real>("mesh.lengthY", 1);
  p.mesh.cellsX = tree.get<int>("mesh.cellsX", 100);
  p.mesh.cellsY = tree.get<int>("mesh.cellsY", 100);

  // Calculate mesh distances
  p.mesh.mesh_dx = p.mesh.sizeX / p.mesh.cellsX;
  p.mesh.mesh_dy = p.mesh.sizeY / p.mesh.cellsY;

  p.sim.final_time = tree.get<Real>("simulation.final_time", 100);
  p.sim.reynolds = tree.get<Real>("simulation.reynolds", 1000);

  return p;
}
