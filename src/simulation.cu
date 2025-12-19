#include "simulation.hpp"

#include "fields/flow_field.hpp"
#include "simulation_impl.cuh"

#include <utility>

Simulation::Simulation(FlowField &field)
    : impl_(new SimulationImpl(field))
{
}

Simulation::~Simulation() { delete impl_; }

void Simulation::run() { impl_->run(); }
