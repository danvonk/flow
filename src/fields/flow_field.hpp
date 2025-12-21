#pragma once

#include "field.hpp"
#include "parameters.hpp"

// Device side view into flow field
struct FlowFieldView {
  // not including boundary cells
  int Nx;
  int Ny;

  IntScalarField::IntScalarFieldView obstacles;
  ScalarField::ScalarFieldView p;
  VectorField::VectorFieldView v;

  VectorField::VectorFieldView fgh;
  ScalarField::ScalarFieldView rhs;
};

class FlowField {
public:
  FlowField(int Nx, int Ny);
  FlowField(const params::Parameters *params);
  virtual ~FlowField() = default;

  inline FlowFieldView view() const
  {
    return {Nx_,
            Ny_,
            obstacles_.view(),
            pressure_.view(),
            velocity_.view(),
            FGH_.view(),
            RHS_.view()};
  };

  // number of field cells (no boundary)
  inline auto Nx() const { return Nx_; }
  inline auto Ny() const { return Ny_; }

  // cells includes the ghost cells here
  inline auto cellsX() const { return cellsX_; }
  inline auto cellsY() const { return cellsY_; }

  inline auto params() const { return params_; }

  inline auto velocity() const { return &velocity_; }
  inline auto pressure() const { return &pressure_; }

protected:
  int Nx_;
  int Ny_;
  // cells include boundaries
  int cellsX_;
  int cellsY_;
  params::Parameters const *params_;

  IntScalarField obstacles_;

  ScalarField pressure_;
  VectorField velocity_;

  VectorField FGH_;
  ScalarField RHS_; // RHS for possion equation
};
