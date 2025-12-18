#include <iostream>

#include "fields/flow_field.hpp"
#include "parameters.hpp"
#include "simulation.hpp"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <cuda_runtime.h>
#include <spdlog/spdlog.h>

namespace po = boost::program_options;

auto main(int argc, char *argv[]) -> int
{

  // command line options and config parsing
  po::options_description desc("Turbulent Flow Simulator Options");
  desc.add_options()("help", "Prints this message of dubious value.")(
      "config", po::value<std::string>(),
      "Path to the .ini config file holding the parameters for the "
      "simulation.");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << '\n';
  }

  params::Parameters params;
  if (vm.count("config")) {
    spdlog::info("Using config file {}", vm["config"].as<std::string>());
    params = params::load_config(vm["config"].as<std::string>());
  }

  // init cuda
  int n_devs = 0;
  cudaGetDeviceCount(&n_devs);
  if (n_devs > 0) {
    cudaDeviceProp p;
    cudaGetDeviceProperties(&p, 0);
    cudaSetDevice(0);

    spdlog::info("Initialised CUDA device 0: {}", p.name);
  }
  else {
    spdlog::error("No CUDA devices found.");
    return -1;
  }

  // create the simulation

  FlowField ff(params);
  Simulation flow_sim(ff);

  flow_sim.initFlowField();

  Real time = 0.0;
  // while (time < params.sim.final_time) {
  flow_sim.solveTimestep();
  time += params.timestep.dt;

  // break;
  // }

  return 0;
}
