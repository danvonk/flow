
#include "fields/flow_field.hpp"
#include "parameters.hpp"
#include "simulation.hpp"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <cuda_runtime.h>
#include <spdlog/spdlog.h>

#include <cfenv>
#include <iostream>


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

  spdlog::info("Flow Simulator. Copyright (c) 2025 Dan Vonk");

  // init cuda
  int n_devs = 0;
  cudaGetDeviceCount(&n_devs);
  if (n_devs > 0) {
    cudaDeviceProp p;
    cudaGetDeviceProperties(&p, 0);
    cudaSetDevice(0);

    spdlog::info("Initialised CUDA device(s):");
    spdlog::info("\t 0: {}", p.name);
  }
  else {
    spdlog::error("No CUDA devices found.");
    return -1;
  }

  params::Parameters params; // host copy
  std::string config_file;
  if (vm.count("config")) {
    spdlog::info("Using config file {}", vm["config"].as<std::string>());
    config_file = vm["config"].as<std::string>();
  }
  params = params::load_config(config_file);

  // create the simulation
  FlowField ff(&params);
  spdlog::info("Created flow field of size {}x{}", ff.Nx(), ff.Ny());
  Simulation flow_sim(ff);

  flow_sim.run();

  return 0;
}
