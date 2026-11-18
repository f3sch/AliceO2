// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DenseSolve.h
/// \brief Small dense factorisations: Cholesky, LU with partial pivoting, and
///        the factor of a positive semi-definite matrix.
///
/// Every routine works in place on a row-major view (data, row stride) of an
/// n x n matrix and allocates nothing: the caller owns all storage. The size
/// is a runtime `int` or a compile-time `std::integral_constant`, so one code
/// serves a caller of fixed size (its loops then of known length) and one
/// whose size is known only at run time. Each loop runs over n: the cost is
/// bounded by n^3 for any input.

#ifndef GTRK_TRACKING_DENSESOLVE_H_
#define GTRK_TRACKING_DENSESOLVE_H_

#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

namespace o2::gtrk::dense
{

/// A row-major matrix in memory the caller owns: element (row, column) at
/// `data[row * stride + column]`.
template <typename T>
struct MatrixView {
  T* data = nullptr;
  int stride = 0;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() T& operator()(int row, int column) const { return data[(row * stride) + column]; }
};

/// Factor the symmetric positive definite \a a, given by its lower triangle,
/// in place: the lower triangle becomes L, a = L L^T. False if a pivot is not
/// positive (the matrix not positive definite, or a NaN); \a a is then
/// partly overwritten.
template <typename T, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool choleskyFactor(MatrixView<T> a, Size size)
{
  const int n = static_cast<int>(size);
  for (int row = 0; row < n; ++row) {
    for (int column = 0; column <= row; ++column) {
      T value = a(row, column);
      for (int k = 0; k < column; ++k) {
        value -= a(row, k) * a(column, k);
      }
      if (row == column) {
        if (!(value > T(0))) {
          return false;
        }
        a(row, row) = sqrtOf(value);
      } else {
        a(row, column) = value / a(column, column);
      }
    }
  }
  return true;
}

/// x <- L^-1 x, L the lower triangle of \a l (forward substitution). With a
/// Cholesky factor this whitens: (L^-1 x)^T (L^-1 y) = x^T a^-1 y.
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void forwardSolve(MatrixView<T> l, Size size, T* x)
{
  const int n = static_cast<int>(size);
  for (int row = 0; row < n; ++row) {
    T value = x[row];
    for (int column = 0; column < row; ++column) {
      value -= l(row, column) * x[column];
    }
    x[row] = value / l(row, row);
  }
}

/// x <- L^-T x, L the lower triangle of \a l (backward substitution).
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void backwardSolve(MatrixView<T> l, Size size, T* x)
{
  const int n = static_cast<int>(size);
  for (int row = n - 1; row >= 0; --row) {
    T value = x[row];
    for (int k = row + 1; k < n; ++k) {
      value -= l(k, row) * x[k];
    }
    x[row] = value / l(row, row);
  }
}

/// x <- a^-1 x, a = L L^T factored by `choleskyFactor` into \a l.
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void choleskySolve(MatrixView<T> l, Size size, T* x)
{
  forwardSolve(l, size, x);
  backwardSolve(l, size, x);
}

/// ln det a, a = L L^T factored by `choleskyFactor` into \a l.
template <typename T, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T choleskyLogDeterminant(MatrixView<T> l, Size size)
{
  const int n = static_cast<int>(size);
  T sum = T(0);
  for (int i = 0; i < n; ++i) {
    sum += logOf(l(i, i));
  }
  return T(2) * sum;
}

/// \a inverse <- a^-1, a = L L^T factored by `choleskyFactor` into \a l;
/// \a inverse must not overlap \a l. Row j of the symmetric inverse is its
/// column j: a^-1 e_j, solved in place in that row.
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void choleskyInverse(MatrixView<T> l, Size size, MatrixView<T> inverse)
{
  const int n = static_cast<int>(size);
  for (int j = 0; j < n; ++j) {
    T* column = &inverse(j, 0);
    for (int i = 0; i < n; ++i) {
      column[i] = i == j ? T(1) : T(0);
    }
    choleskySolve(l, size, column);
  }
}

/// Factor \a a in place as P a = L U with partial pivoting (L unit lower,
/// below the diagonal; U on and above it): column k's pivot is its entry of
/// largest magnitude on or below the diagonal, rows k and \a pivots[k]
/// swapped. False if a pivot is exactly zero (the matrix singular).
template <typename T, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool luFactor(MatrixView<T> a, Size size, int* pivots)
{
  const int n = static_cast<int>(size);
  for (int k = 0; k < n; ++k) {
    int pivot = k;
    T largest = o2::gpu::GPUCommonMath::Abs(a(k, k));
    for (int row = k + 1; row < n; ++row) {
      const T magnitude = o2::gpu::GPUCommonMath::Abs(a(row, k));
      if (magnitude > largest) {
        largest = magnitude;
        pivot = row;
      }
    }
    pivots[k] = pivot;
    if (!(largest != T(0))) {
      return false;
    }
    if (pivot != k) {
      for (int column = 0; column < n; ++column) {
        const T swapped = a(k, column);
        a(k, column) = a(pivot, column);
        a(pivot, column) = swapped;
      }
    }
    for (int row = k + 1; row < n; ++row) {
      const T factor = a(row, k) / a(k, k);
      a(row, k) = factor;
      for (int column = k + 1; column < n; ++column) {
        a(row, column) -= factor * a(k, column);
      }
    }
  }
  return true;
}

/// x <- a^-1 x, a factored by `luFactor` into \a lu and \a pivots.
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void luSolve(MatrixView<T> lu, Size size, const int* pivots, T* x)
{
  const int n = static_cast<int>(size);
  for (int k = 0; k < n; ++k) {
    if (pivots[k] != k) {
      const T swapped = x[k];
      x[k] = x[pivots[k]];
      x[pivots[k]] = swapped;
    }
  }
  for (int row = 1; row < n; ++row) {
    T value = x[row];
    for (int column = 0; column < row; ++column) {
      value -= lu(row, column) * x[column];
    }
    x[row] = value;
  }
  for (int row = n - 1; row >= 0; --row) {
    T value = x[row];
    for (int column = row + 1; column < n; ++column) {
      value -= lu(row, column) * x[column];
    }
    x[row] = value / lu(row, row);
  }
}

/// ln |det a|, a factored by `luFactor` into \a lu.
template <typename T, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T luLogAbsDeterminant(MatrixView<T> lu, Size size)
{
  const int n = static_cast<int>(size);
  T sum = T(0);
  for (int i = 0; i < n; ++i) {
    sum += logOf(o2::gpu::GPUCommonMath::Abs(lu(i, i)));
  }
  return sum;
}

/// A factor F of the positive semi-definite \a a (n x n, both triangles
/// given), F F^T = a, into \a factor: the symmetrically pivoted L D L^T,
/// F = P^T L D^1/2. Each step pivots on the largest diagonal of what is left;
/// once that is not positive, what is left is zero (a semi-definite matrix
/// with a zero diagonal is zero; below it only by rounding) and F's remaining
/// columns are zero. F is not triangular. \a a is overwritten; \a order holds
/// n ints of scratch.
template <typename T, typename Size>
GPUhdi() GPUAlwaysInline() void semidefiniteFactor(MatrixView<T> a, Size size, int* order, MatrixView<T> factor)
{
  const int n = static_cast<int>(size);
  for (int i = 0; i < n; ++i) {
    order[i] = i;
  }
  int rank = 0;
  for (int k = 0; k < n; ++k) {
    int pivot = k;
    for (int i = k + 1; i < n; ++i) {
      if (a(i, i) > a(pivot, pivot)) {
        pivot = i;
      }
    }
    if (!(a(pivot, pivot) > T(0))) {
      break;
    }
    if (pivot != k) {
      for (int j = 0; j < n; ++j) {
        const T swapped = a(k, j);
        a(k, j) = a(pivot, j);
        a(pivot, j) = swapped;
      }
      for (int i = 0; i < n; ++i) {
        const T swapped = a(i, k);
        a(i, k) = a(i, pivot);
        a(i, pivot) = swapped;
      }
      const int swapped = order[k];
      order[k] = order[pivot];
      order[pivot] = swapped;
    }
    const T d = a(k, k);
    for (int i = k + 1; i < n; ++i) {
      a(i, k) /= d;
    }
    for (int i = k + 1; i < n; ++i) {
      for (int j = k + 1; j < n; ++j) {
        a(i, j) -= a(i, k) * a(j, k) * d;
      }
    }
    rank = k + 1;
  }
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      T value = T(0);
      if (k < rank && i >= k) {
        value = (i == k ? T(1) : a(i, k)) * sqrtOf(a(k, k));
      }
      factor(order[i], k) = value;
    }
  }
}

} // namespace o2::gtrk::dense

#endif // GTRK_TRACKING_DENSESOLVE_H_
