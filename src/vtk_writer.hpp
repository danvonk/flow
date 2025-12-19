#pragma once

#include "config.hpp"
#include "fields/flow_field.hpp"

#include <vtk/vtkHDFReader.h>
#include <vtk/vtkIOHDFModule.h>
#include <vtk/vtkImageData.h>

class VTKWriter {
public:
  VTKWriter();

  void write_flow_field(FlowField &field);

private:
  vtkHDFWriter wr_;
};
