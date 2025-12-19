add_executable(tflow)

target_sources(tflow
  PRIVATE
  simulation.hpp
  simulation.cu
  simulation_impl.cuh
  simulation_impl.cu
  parameters.hpp
  parameters.cu
  vtk_writer.hpp
  vtk_writer.cpp
  config.hpp
  fields/flow_field.hpp
  fields/flow_field.cpp
  fields/field.hpp
  stencils/stencil.hpp
  stencils/iterator.cuh
  stencils/rhs_stencil.cuh
  stencils/fgh_stencil.cuh
  stencils/moving_wall_stencil.cuh
  main.cpp)


target_link_libraries(tflow PRIVATE
  spdlog::spdlog
  CUDA::cudart
  Boost::program_options
  VTK::IOHDF
  VTK::IOXML
  VTK::IOCore
  VTK::CommonDataModel
  VTK::CommonCore
  Boost::boost # target for header-only libraries
)

set_target_properties(tflow PROPERTIES
CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(tflow
PRIVATE
${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS tflow)

