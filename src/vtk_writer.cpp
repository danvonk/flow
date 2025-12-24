#include "vtk_writer.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>

void VTKWriter::write_flow_field(FlowField &field, Real time)
{
  std::vector<Real> velocity_data;
  std::vector<Real> pressure_data;

  field.velocity()->to_host(velocity_data);
  field.pressure()->to_host(pressure_data);

  const int NxTot = field.cellsX();
  const int NyTot = field.cellsY();

  // for (int j = 2; j < NyTot - 1; ++j) {
  //   for (int i = 2; i < NxTot - 1; ++i) {
  //     auto here = 2 * (i + NxTot * j);
  //     auto left = 2 * ((i - 1) + NxTot * j);
  //     auto down = 2 * (i + NxTot * (j - 1));
  //     std::cout << velocity_data[here + 0] << "," << velocity_data[here + 1]
  //               << " ";
  //   }
  //   std::cout << '\n';
  // }

  vtkNew<vtkImageData> image;
  image->SetOrigin(0, 0, 0);
  image->SetSpacing(field.params()->mesh.mesh_dx, field.params()->mesh.mesh_dy,
                    1.0);
  image->SetDimensions(field.Nx() + 1, field.Ny() + 1, 2);

  vtkNew<vtkDoubleArray> velocity;
  velocity->SetName("velocity");
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(field.Nx() * field.Ny());

  auto *ptr = velocity->GetPointer(0);
  for (int j = 2; j < NyTot - 1; ++j) {
    for (int i = 2; i < NxTot - 1; ++i) {
      auto here = 2 * (i + NxTot * j);
      auto left = 2 * ((i - 1) + NxTot * j);
      auto down = 2 * (i + NxTot * (j - 1));

      auto u = 0.5 * (velocity_data[here + 0] + velocity_data[left + 0]);
      auto v = 0.5 * (velocity_data[here + 1] + velocity_data[down + 1]);

      const auto cell = (i - 2) + field.Nx() * (j - 2);

      ptr[3 * cell + 0] = u;
      ptr[3 * cell + 1] = v;
      ptr[3 * cell + 2] = 0.0;
    }
  }

  image->GetCellData()->AddArray(velocity);

  vtkNew<vtkDoubleArray> pressure;
  pressure->SetName("pressure");
  pressure->SetNumberOfComponents(1);
  pressure->SetNumberOfTuples(field.Nx() * field.Ny());

  ptr = pressure->GetPointer(0);
  for (int j = 2; j < NyTot - 1; ++j) {
    for (int i = 2; i < NxTot - 1; ++i) {
      auto here = 2 * (i + NxTot * j);
      const auto cell = (i - 2) + field.Nx() * (j - 2);
      ptr[cell] = here;
    }
  }
  image->GetCellData()->AddArray(pressure);

  // TODO: use std::format etc.
  char buffer[256];
  sprintf(buffer, "flow%.0f.vti", time * 100);

  // int dims[3];
  // image->GetDimensions(dims);
  // auto nCells = image->GetNumberOfCells();
  // std::cout << "dims: " << dims[0] << "," << dims[1] << "," << dims[2] <<
  // "\n"; std::cout << "cells: " << nCells << " tuples: " <<
  // velocity->GetNumberOfTuples() << "\n";

  writer_->SetFileName(buffer);
  writer_->SetInputData(image);
  writer_->Write();
}
