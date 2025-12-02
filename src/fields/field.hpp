#pragma once

#include <cuda.h>
#include <cuda_runtime.h>

// Fields are stored directly on GPU memory
template <typename T = float> class Field {
public:
  Field(int Nx, int Ny, int components)
      : sizeX_(Nx),
        sizeY_(Ny),
        components_(components), // scalar, 2D, 3D vector etc.
        size_(components * Nx * Ny * sizeof(T))
  {
    cudaMalloc(&data_, size_);
  };

  virtual ~Field() { cudaFree(data_); }

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

class ScalarField : public Field<float> {
public:
  ScalarField(int Nx, int Ny)
      : Field(Nx, Ny, 1)
  {
  }

  float &getScalar(int i, int j);
};

class IntScalarField : public Field<int> {
public:
  IntScalarField(int Nx, int Ny)
      : Field(Nx, Ny, 1)
  {
  }

  float &getScalar(int i, int j);
};

class VectorField : public Field<> {
public:
  VectorField(int Nx, int Ny, int components)
      : Field(Nx, Ny, components)
  {
  }

  float *getVector(int i, int j);
};
