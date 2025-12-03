#include "parameters.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

using params::Parameters;

void Parameters::load_config(const std::string &file_path)
{
  if (file_path.empty())
    return;

  pt::ptree tree;
  pt::read_ini(file_path, tree);

  timestep.dt = tree.get<Real>("timestep.dt", 1.0);
  timestep.tau = tree.get<Real>("timestep.tau", 1.0);

  mesh.sizeX = tree.get<Real>("mesh.lengthX", 1);
  mesh.sizeY = tree.get<Real>("mesh.lengthY", 1);
  mesh.cellsX = tree.get<Real>("mesh.cellsX", 100);
  mesh.cellsY = tree.get<Real>("mesh.cellsY", 100);

  sim.final_time = tree.get<Real>("simulation.final_time", 100);
  sim.reynolds = tree.get<Real>("simulation.reynolds", 1000);
}
