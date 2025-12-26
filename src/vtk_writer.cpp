#include "vtk_writer.hpp"

#include <cassert>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHDFWriter.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>


#include <spdlog/spdlog.h>

void VTKWriter::write_flow_field(FlowField& field, int step)
{
  std::vector<Real> velocity_data;
  std::vector<Real> pressure_data;

  field.velocity()->to_host(velocity_data);
  field.pressure()->to_host(pressure_data);

  const int NxTot = field.cellsX();
  const int NyTot = field.cellsY();

  assert(pressure_data.size() == size_t(NxTot * NyTot));
  assert(velocity_data.size() == size_t(2 * NxTot * NyTot));

  const int nx = field.Nx();
  const int ny = field.Ny();

  const double dx = field.params()->mesh.mesh_dx;
  const double dy = field.params()->mesh.mesh_dy;

  vtkNew<vtkPolyData> poly;
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints((nx + 1) * (ny + 1));

  auto pid = [nx](int i, int j) -> vtkIdType {
    return static_cast<vtkIdType>(i + (nx + 1) * j);
  };

  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i <= nx; ++i) {
      points->SetPoint(pid(i, j), i * dx, j * dy, 0.0);
    }
  }

  vtkNew<vtkCellArray> polys;
  polys->AllocateExact(nx * ny, nx * ny * 4);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      vtkIdType ids[4] = {
        pid(i,     j),
        pid(i + 1, j),
        pid(i + 1, j + 1),
        pid(i,     j + 1)
      };
      polys->InsertNextCell(4, ids);
    }
  }

  poly->SetPoints(points);
  poly->SetPolys(polys);

  const auto nCells = poly->GetNumberOfCells();
  assert(nCells == (nx * ny));
                        
  vtkNew<vtkDoubleArray> velocity;
  velocity->SetName("velocity");
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(nCells);

  vtkNew<vtkDoubleArray> pressure;
  pressure->SetName("pressure");
  pressure->SetNumberOfComponents(1);
  pressure->SetNumberOfTuples(nCells);

  auto* vptr = velocity->GetPointer(0);
  auto* pptr = pressure->GetPointer(0);

  for (int j = 2; j < NyTot - 1; ++j) {
    for (int i = 2; i < NxTot - 1; ++i) {

      const int cell = (i - 2) + nx * (j - 2);
      assert(cell >= 0 && cell < nx * ny);

      const int here_v = 2 * (i + NxTot * j);
      const int left_v = 2 * ((i - 1) + NxTot * j);
      const int down_v = 2 * (i + NxTot * (j - 1));

      const double u = 0.5 * (velocity_data[here_v + 0] + velocity_data[left_v + 0]);
      const double v = 0.5 * (velocity_data[here_v + 1] + velocity_data[down_v + 1]);

      vptr[3 * cell + 0] = u;
      vptr[3 * cell + 1] = v;
      vptr[3 * cell + 2] = 0.0;

      const int here_p = i + NxTot * j;
      pptr[cell + 0] = pressure_data[here_p];

    }
  }

  auto* cd = poly->GetCellData();
  cd->AddArray(velocity);
  cd->AddArray(pressure);
  cd->SetActiveVectors("velocity");
  cd->SetActiveScalars("pressure");

  char buffer[256];
  sprintf(buffer, "flow_%06d.vtkhdf", step);

  vtkNew<vtkHDFWriter> writer;
  writer->SetFileName(buffer);
  writer->SetOverwrite(true);
  writer->SetInputData(poly);
  writer->Write();
}

