#pragma once

#include "config.hpp"
#include "fields/flow_field.hpp"

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

class VTKWriter {
public:
  VTKWriter() {};

  void write_flow_field(FlowField &field, Real dt);

private:
  vtkNew<vtkXMLImageDataWriter> writer_;
};
