#include "flow_field.hpp"
#include "fields/field.hpp"
#include "parameters.hpp"

FlowField::FlowField(int Nx, int Ny)
    : Nx_(Nx),
      Ny_(Ny),
      // extra 2 ghost cells on left and 1 on right
      cellsX_(Nx + 3),
      // extra 2 ghost cells on bottom and 1 on top
      cellsY_(Ny + 3),
      obstacles_(IntScalarField(cellsX_, cellsY_)),
      pressure_(ScalarField(cellsX_, cellsY_)),
      velocity_(VectorField(cellsX_, cellsY_)),
      FGH_(VectorField(cellsX_, cellsY_)),
      RHS_(ScalarField(cellsX_, cellsY_))
{
}

FlowField::FlowField(const params::Parameters *params)
    : Nx_(params->mesh.cellsX),
      Ny_(params->mesh.cellsY),
      // extra 2 ghost cells on left and 1 on right
      cellsX_(Nx_ + 3),
      // extra 2 ghost cells on bottom and 1 on top
      cellsY_(Ny_ + 3),
      params_(params),
      obstacles_(IntScalarField(cellsX_, cellsY_)),
      pressure_(ScalarField(cellsX_, cellsY_)),
      velocity_(VectorField(cellsX_, cellsY_)),
      FGH_(VectorField(cellsX_, cellsY_)),
      RHS_(ScalarField(cellsX_, cellsY_))
{
}
