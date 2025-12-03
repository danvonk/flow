#pragma once

#include "config.hpp"

namespace stencils {
class RHSStencil {
public:
  __device__ void operator()(FlowField &field, int i, int j) const { return; }
}
} // namespace stencils
