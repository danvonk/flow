#include "simulation.hpp"

#include "fields/flow_field.hpp"
#include "simulation_impl.cuh"

#include <utility>

Simulation::Simulation(FlowField &field)
    : impl_(new SimulationImpl(field))
{
}

Simulation::~Simulation() { delete impl_; }

void Simulation::initFlowField() { impl_->initFlowField(); }

void Simulation::solveTimestep() { impl_->solveTimestep(); }

void Simulation::setTimestep() { impl_->setTimestep(); }
