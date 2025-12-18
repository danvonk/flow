add_executable(tflow)

target_sources(tflow
  PRIVATE
  simulation.hpp
  simulation.cu
  simulation_impl.cuh
  simulation_impl.cu
  parameters.hpp
  parameters.cpp
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
  hdf5::hdf5_cpp
  hdf5::hdf5
  hdf5::hdf5_hl_cpp
  hdf5::hdf5_hl
  Boost::program_options
  Boost::boost # target for header-only libraries
)

target_include_directories(tflow
PRIVATE
${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS tflow)

# target_sources(simple
#     PUBLIC
#     $<INSTALL_INTERFACE:include/src/simple.h>
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/simple.h>
#     PRIVATE simple.cc)

# if(ENABLE_TESTING)
#     MAKE_CC_TEST(
#         NAME simple_test
#         SRCS simple_test.cc
#         DEPS simple)
# endif()
