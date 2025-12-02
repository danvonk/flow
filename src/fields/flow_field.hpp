#pragma once

#include "field.hpp"

class FlowField {
public:
  FlowField(int Nx, int Ny);
  virtual ~FlowField() = default;

protected:
  ScalarField pressure_;
  VectorField velocity_;

  VectorField FGH_;
  ScalarField RHS_; // RHS for possion equation
};
