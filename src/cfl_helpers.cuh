#pragma once

#include <fields/flow_field.hpp>

Real calc_new_timestep(FlowFieldView field, cudaStream_t stream = 0);
