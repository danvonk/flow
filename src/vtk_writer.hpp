#pragma once

#include "config.hpp"
#include "fields/flow_field.hpp"

#include <vtkHDFWriter.h>
#include <vtkNew.h>

class VTKWriter {
public:
  VTKWriter() {};

  void write_flow_field(FlowField &field, int step);

private:
  vtkNew<vtkHDFWriter> writer_;
};
