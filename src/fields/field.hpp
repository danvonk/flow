#pragma once

#include "config.hpp"

#include <driver_types.h>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

// Fields are stored directly on GPU memory but their memory is managed on the
// CPU through normal C++ RAII
template <typename T = Real> class Field {
public:
  Field(int Nx, int Ny, int components)
      : sizeX_(Nx),
        sizeY_(Ny),
        components_(components), // scalar, 2D, 3D vector etc.
        size_(components * Nx * Ny * sizeof(T))
  {
    cudaMalloc(&data_, size_);
    cudaMemset(data_, 0.0, size_);
  };

  virtual ~Field() { cudaFree(data_); }

  void to_host(std::vector<T> &host) const
  {
    host.resize(sizeX_ * sizeY_ * components_);
    cudaMemcpy(host.data(), data_, host.size(), cudaMemcpyDeviceToHost);
  }

  Field(const Field &) = delete;
  Field &operator=(const Field &) = delete;

  auto getNx() const { return sizeX_; }
  auto getNy() const { return sizeY_; }

protected:
  T *data_ = nullptr;
  int sizeX_;
  int sizeY_;
  int components_;
  int size_;
};

class ScalarField : public Field<Real> {
public:
  ScalarField(int Nx, int Ny)
      : Field(Nx, Ny, 1)
  {
  }

  // GPU-side view of the field
  struct ScalarFieldView {
    Real *data;
    int Nx;
    int Ny;

    HD inline Real &get(int i, int j) { return data[i + Nx * j]; }
  };

  ScalarFieldView view() const { return {data_, sizeX_, sizeY_}; };
};

class IntScalarField : public Field<int> {
public:
  IntScalarField(int Nx, int Ny)
      : Field(Nx, Ny, 1)
  {
  }

  struct IntScalarFieldView {
    int *data;
    int Nx;
    int Ny;
  };

  inline IntScalarFieldView view() const { return {data_, sizeX_, sizeY_}; };
};

class VectorField : public Field<> {
public:
  VectorField(int Nx, int Ny)
      : Field(Nx, Ny, 2)
  {
  }

  struct VectorFieldView {
    Real *data;
    int Nx;
    int Ny;

    HD inline Real &u(int i, int j) { return data[2 * (i + Nx * j) + 0]; }

    HD inline Real &v(int i, int j) { return data[2 * (i + Nx * j) + 1]; }
  };

  VectorFieldView view() const { return {data_, sizeX_, sizeY_}; };
};
