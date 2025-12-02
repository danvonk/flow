add_executable(tflow)

target_sources(tflow
  PRIVATE
  simulation.hpp
  parameters.hpp
  config.hpp
  fields/flow_field.hpp
  fields/field.hpp
  stencils/rhs_stencil.cuh
  main.cpp)


target_link_libraries(tflow PRIVATE
  spdlog::spdlog
  CUDA::cudart
  CUDA::cufft
  hdf5::hdf5_cpp
  hdf5::hdf5
  hdf5::hdf5_hl_cpp
  hdf5::hdf5_hl
)

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
