#include "vtk_writer.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHDFWriter.h>
#include <vtkImageData.h>

void VTKWriter::write_flow_field(FlowField &field, Real dt)
{
  std::vector<Real> velocity_data;
  std::vector<Real> pressure_data;

  field.velocity()->to_host(velocity_data);
  field.pressure()->to_host(pressure_data);

  for (auto j = 0; j < field.Ny(); ++j) {
    for (auto i = 0; i < field.Nx(); ++i) {
      std::cout << velocity_data[2 * (i + (field.Nx() + 3) * j) + 0] << ","
                << velocity_data[2 * (i + (field.Nx() + 3) * j) + 1] << " ";
    }
    std::cout << '\n';
  }

  vtkNew<vtkImageData> image;
  image->SetOrigin(0, 0, 0);
  image->SetSpacing(field.params()->mesh.mesh_dx, field.params()->mesh.mesh_dy,
                    1.0);
  image->SetDimensions(field.Nx() + 1, field.Ny() + 1, 1);

  vtkNew<vtkDoubleArray> velocity;
  velocity->SetName("velocity");
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(field.Nx() * field.Ny());

  const int NxTot = field.Nx() + 3;
  const int NyTot = field.Ny() + 3;

  auto *ptr = velocity->GetPointer(0);
  for (int j = 2; j <= field.Ny() + 1; ++j) {
    for (int i = 2; i <= field.Nx() + 1; ++i) {
      auto here = 2 * (i + NxTot * j);
      auto left = 2 * ((i - 1) + NxTot * j);
      auto down = 2 * (i + NxTot * (j - 1));

      auto u = 0.5 * (velocity_data[here + 0] + velocity_data[left + 0]);
      auto v = 0.5 * (velocity_data[here + 1] + velocity_data[down + 1]);

      const auto cell = (i - 2) + field.Nx() * (j - 2);

      ptr[3 * cell + 0] = u;
      ptr[3 * cell + 1] = v;
      ptr[3 * cell + 2] = 0.;
    }
  }

  image->GetCellData()->AddArray(velocity);

  char buffer[256];
  sprintf(buffer, "flow%f.vti", dt);

  writer_->SetFileName(buffer);
  writer_->SetInputData(image);
  writer_->Write();
}
