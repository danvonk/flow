#pragma once

#include "field.hpp"
#include "parameters.hpp"

// Device side view into flow field
struct FlowFieldView {
  int cellsX;
  int cellsY;

  IntScalarField::IntScalarFieldView obstacles;
  ScalarField::ScalarFieldView p;
  VectorField::VectorFieldView v;

  VectorField::VectorFieldView fgh;
  ScalarField::ScalarFieldView rhs;

  params::Parameters params;
};

class FlowField {
public:
  FlowField(int Nx, int Ny);
  FlowField(const params::Parameters &params);
  virtual ~FlowField() = default;

  inline FlowFieldView view() const
  {
    return {cellsX_,          cellsY_,     obstacles_.view(), pressure_.view(),
            velocity_.view(), FGH_.view(), RHS_.view(),       params_};
  };

  inline auto cellsX() const { return cellsX_; }

  inline auto cellsY() const { return cellsY_; }

protected:
  int Nx_;
  int Ny_;
  // cells include boundaries
  int cellsX_;
  int cellsY_;
  params::Parameters params_;

  IntScalarField obstacles_;

  ScalarField pressure_;
  VectorField velocity_;

  VectorField FGH_;
  ScalarField RHS_; // RHS for possion equation
};
