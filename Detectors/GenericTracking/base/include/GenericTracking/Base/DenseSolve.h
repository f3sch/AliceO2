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
/// \brief Small dense factorisations of symmetric matrices: Cholesky, the
///        symmetric indefinite L D L^T (Bunch-Kaufman), and the factor of a
///        positive semi-definite matrix.
///
/// Every routine works in place on a view of an n x n matrix in memory the
/// caller owns and allocates nothing. A view is a row-major `MatrixView` or a
/// `PackedView`, a symmetric matrix by its lower triangle; every routine
/// reads and writes the lower triangle only, so both serve. The size is a
/// runtime `int` or a compile-time `std::integral_constant`, so one code
/// serves a caller of fixed size (its loops then of known length) and one
/// whose size is known only at run time. Each loop runs over n: the cost is
/// bounded by n^3 for any input.

#ifndef GTRK_BASE_DENSESOLVE_H_
#define GTRK_BASE_DENSESOLVE_H_

#include "GenericTracking/Base/ElementaryMath.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

#include <type_traits>

namespace o2::gtrk::dense
{

/// A row-major matrix: element (row, column) at `data[row * stride + column]`.
template <typename T>
struct MatrixView {
  T* data = nullptr;
  int stride = 0;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() T& operator()(int row, int column) const { return data[(row * stride) + column]; }
};

/// Where (row, column) of a symmetric matrix sits in its packed lower
/// triangle: row by row, row r's columns 0..r at r (r + 1) / 2.
[[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr int packedIndex(int row, int column)
{
  return row >= column ? ((row * (row + 1)) / 2) + column : ((column * (column + 1)) / 2) + row;
}

/// The entries a packed symmetric n x n matrix holds, n (n + 1) / 2.
[[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr int packedSize(int n) { return (n * (n + 1)) / 2; }

/// A symmetric matrix by its lower triangle, packed (`packedIndex`):
/// (row, column) and (column, row) are the same element.
template <typename T>
struct PackedView {
  T* data = nullptr;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() T& operator()(int row, int column) const { return data[packedIndex(row, column)]; }
};

/// A view's element type.
template <typename View>
using ElementOf = std::remove_reference_t<decltype(std::declval<View>()(0, 0))>;

/// Factor the symmetric positive definite \a a in place: its lower triangle
/// becomes L, a = L L^T. False if a pivot is not positive (the matrix not
/// positive definite, or a NaN); \a a is then partly overwritten.
template <typename View, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool choleskyFactor(View a, Size size)
{
  using T = ElementOf<View>;
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
template <typename View, typename Size, typename T>
GPUhdi() GPUAlwaysInline() void forwardSolve(View l, Size size, T* x)
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
template <typename View, typename Size, typename T>
GPUhdi() GPUAlwaysInline() void backwardSolve(View l, Size size, T* x)
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
template <typename View, typename Size, typename T>
GPUhdi() GPUAlwaysInline() void choleskySolve(View l, Size size, T* x)
{
  forwardSolve(l, size, x);
  backwardSolve(l, size, x);
}

/// ln det a, a = L L^T factored by `choleskyFactor` into \a l.
template <typename View, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() ElementOf<View> choleskyLogDeterminant(View l, Size size)
{
  using T = ElementOf<View>;
  const int n = static_cast<int>(size);
  T sum = T(0);
  for (int i = 0; i < n; ++i) {
    sum += logOf(l(i, i));
  }
  return T(2) * sum;
}

/// \a inverse <- a^-1, a = L L^T factored by `choleskyFactor` into \a l;
/// \a inverse must not overlap \a l. First M = L^-1 into the lower triangle,
/// then a^-1 = M^T M in place, (i, j) = sum over k >= i of M(k, i) M(k, j)
/// for i >= j, column by column and down each column: an entry overwritten
/// is one no later entry reads. A row-major \a inverse gets its upper
/// triangle mirrored.
template <typename View, typename Size, typename OutView>
GPUhdi() GPUAlwaysInline() void choleskyInverse(View l, Size size, OutView inverse)
{
  using T = ElementOf<OutView>;
  const int n = static_cast<int>(size);
  for (int column = 0; column < n; ++column) {
    inverse(column, column) = T(1) / l(column, column);
    for (int row = column + 1; row < n; ++row) {
      T value = T(0);
      for (int k = column; k < row; ++k) {
        value -= l(row, k) * inverse(k, column);
      }
      inverse(row, column) = value / l(row, row);
    }
  }
  for (int column = 0; column < n; ++column) {
    for (int row = column; row < n; ++row) {
      T value = T(0);
      for (int k = row; k < n; ++k) {
        value += inverse(k, row) * inverse(k, column);
      }
      inverse(row, column) = value;
    }
  }
  if constexpr (std::is_same_v<OutView, MatrixView<T>>) {
    for (int row = 0; row < n; ++row) {
      for (int column = row + 1; column < n; ++column) {
        inverse(row, column) = inverse(column, row);
      }
    }
  }
}

namespace detail
{
/// Exchange rows and columns \a k and \a p (k < p) of the symmetric matrix
/// held by \a a's lower triangle, the columns before k included (the rows
/// of L found so far): the lower triangle of P a P^T.
template <typename View>
GPUhdi() GPUAlwaysInline() void exchange(View a, int n, int k, int p)
{
  using T = ElementOf<View>;
  const auto swap = [](T& x, T& y) {
    const T t = x;
    x = y;
    y = t;
  };
  for (int j = 0; j < k; ++j) {
    swap(a(k, j), a(p, j));
  }
  swap(a(k, k), a(p, p));
  for (int i = k + 1; i < p; ++i) {
    T& ik = a(i, k);
    T& pi = a(p, i);
    swap(ik, pi);
  }
  for (int i = p + 1; i < n; ++i) {
    swap(a(i, k), a(i, p));
  }
}
} // namespace detail

/// Factor the symmetric (possibly indefinite) \a a in place as
/// P a P^T = L D L^T, Bunch and Kaufman's partial pivoting: D block diagonal
/// with 1 x 1 and 2 x 2 blocks, L unit lower. Column k's pivot is its
/// diagonal if that is at least alpha times the column's largest entry below
/// it (or than that relative to the largest entry of the row it sits in);
/// else the row's own diagonal, or the 2 x 2 block of k and that row. alpha =
/// (1 + sqrt 17) / 8 bounds the growth of the entries (Bunch and Kaufman,
/// 1977). The lower triangle holds L below each block and D on it (a 2 x 2
/// block's off-diagonal included); \a swaps[k] is the row exchanged at step
/// k, \a blocks[k] the size of the block starting at k (0 for a 2 x 2
/// block's second row). False if a pivot block is exactly singular.
template <typename View, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool ldltFactor(View a, Size size, int* swaps, int* blocks)
{
  using T = ElementOf<View>;
  const int n = static_cast<int>(size);
  const T alpha = (T(1) + sqrtOf(T(17))) / T(8);
  for (int k = 0; k < n; ++k) {
    swaps[k] = k;
    blocks[k] = 1;
  }
  int k = 0;
  while (k < n) {
    const T diagonal = o2::gpu::GPUCommonMath::Abs(a(k, k));
    int largestRow = k;
    T columnLargest = T(0);
    for (int i = k + 1; i < n; ++i) {
      const T magnitude = o2::gpu::GPUCommonMath::Abs(a(i, k));
      if (magnitude > columnLargest) {
        columnLargest = magnitude;
        largestRow = i;
      }
    }
    if (!(diagonal != T(0) || columnLargest != T(0))) {
      return false;
    }
    int step = 1;
    int pivot = k;
    if (!(diagonal >= alpha * columnLargest)) {
      T rowLargest = T(0);
      for (int j = k; j < n; ++j) {
        if (j != largestRow) {
          rowLargest = o2::gpu::GPUCommonMath::Max(rowLargest, o2::gpu::GPUCommonMath::Abs(j < largestRow ? a(largestRow, j) : a(j, largestRow)));
        }
      }
      if (diagonal * rowLargest >= alpha * columnLargest * columnLargest) {
        pivot = k;
      } else if (o2::gpu::GPUCommonMath::Abs(a(largestRow, largestRow)) >= alpha * rowLargest) {
        pivot = largestRow;
      } else {
        pivot = largestRow;
        step = 2;
      }
    }
    const int swapped = k + step - 1;
    if (pivot != swapped) {
      detail::exchange(a, n, swapped, pivot);
    }
    swaps[swapped] = pivot;
    if (step == 1) {
      const T d = a(k, k);
      for (int i = k + 1; i < n; ++i) {
        const T factor = a(i, k) / d;
        for (int j = k + 1; j <= i; ++j) {
          a(i, j) -= factor * a(j, k);
        }
      }
      for (int i = k + 1; i < n; ++i) {
        a(i, k) /= d;
      }
    } else {
      const T d11 = a(k, k);
      const T d21 = a(k + 1, k);
      const T d22 = a(k + 1, k + 1);
      const T determinant = (d11 * d22) - (d21 * d21);
      if (!(determinant != T(0))) {
        return false;
      }
      for (int i = k + 2; i < n; ++i) {
        // (l_i1, l_i2) = (a_ik, a_ik+1) D^-1
        const T first = ((a(i, k) * d22) - (a(i, k + 1) * d21)) / determinant;
        const T second = ((a(i, k + 1) * d11) - (a(i, k) * d21)) / determinant;
        for (int j = k + 2; j <= i; ++j) {
          a(i, j) -= (first * a(j, k)) + (second * a(j, k + 1));
        }
      }
      for (int i = k + 2; i < n; ++i) {
        const T first = ((a(i, k) * d22) - (a(i, k + 1) * d21)) / determinant;
        const T second = ((a(i, k + 1) * d11) - (a(i, k) * d21)) / determinant;
        a(i, k) = first;
        a(i, k + 1) = second;
      }
      blocks[k] = 2;
      blocks[k + 1] = 0;
    }
    k += step;
  }
  return true;
}

/// x <- a^-1 x, a factored by `ldltFactor` into \a ldl, \a swaps, \a blocks.
template <typename View, typename Size, typename T>
GPUhdi() GPUAlwaysInline() void ldltSolve(View ldl, Size size, const int* swaps, const int* blocks, T* x)
{
  const int n = static_cast<int>(size);
  for (int k = 0; k < n; ++k) {
    if (swaps[k] != k) {
      const T t = x[k];
      x[k] = x[swaps[k]];
      x[swaps[k]] = t;
    }
  }
  // L y = P x
  for (int k = 0; k < n; ++k) {
    const int first = blocks[k] == 0 ? k + 1 : k + blocks[k];
    for (int i = first; i < n; ++i) {
      x[i] -= ldl(i, k) * x[k];
    }
  }
  // D z = y
  for (int k = 0; k < n;) {
    if (blocks[k] == 1) {
      x[k] /= ldl(k, k);
      k += 1;
    } else {
      const T d11 = ldl(k, k);
      const T d21 = ldl(k + 1, k);
      const T d22 = ldl(k + 1, k + 1);
      const T determinant = (d11 * d22) - (d21 * d21);
      const T first = ((d22 * x[k]) - (d21 * x[k + 1])) / determinant;
      const T second = ((d11 * x[k + 1]) - (d21 * x[k])) / determinant;
      x[k] = first;
      x[k + 1] = second;
      k += 2;
    }
  }
  // L^T w = z
  for (int k = n - 1; k >= 0; --k) {
    const int first = blocks[k] == 0 ? k + 1 : k + blocks[k];
    T value = x[k];
    for (int i = first; i < n; ++i) {
      value -= ldl(i, k) * x[i];
    }
    x[k] = value;
  }
  for (int k = n - 1; k >= 0; --k) {
    if (swaps[k] != k) {
      const T t = x[k];
      x[k] = x[swaps[k]];
      x[swaps[k]] = t;
    }
  }
}

/// ln |det a|, a factored by `ldltFactor`: the blocks' determinants.
template <typename View, typename Size>
[[nodiscard]] GPUhdi() GPUAlwaysInline() ElementOf<View> ldltLogAbsDeterminant(View ldl, Size size, const int* blocks)
{
  using T = ElementOf<View>;
  const int n = static_cast<int>(size);
  T sum = T(0);
  for (int k = 0; k < n; ++k) {
    if (blocks[k] == 1) {
      sum += logOf(o2::gpu::GPUCommonMath::Abs(ldl(k, k)));
    } else if (blocks[k] == 2) {
      const T determinant = (ldl(k, k) * ldl(k + 1, k + 1)) - (ldl(k + 1, k) * ldl(k + 1, k));
      sum += logOf(o2::gpu::GPUCommonMath::Abs(determinant));
    }
  }
  return sum;
}

/// A factor F of the positive semi-definite \a a, F F^T = a, into the
/// row-major \a factor: the symmetrically pivoted L D L^T, F = P^T L D^1/2.
/// Each step pivots on the largest diagonal of what is left; once that is not
/// positive, what is left is zero (a semi-definite matrix with a zero
/// diagonal is zero; below it only by rounding) and F's remaining columns are
/// zero. F is not triangular. \a a (its lower triangle) is overwritten;
/// \a order holds n ints of scratch.
template <typename View, typename Size, typename T>
GPUhdi() GPUAlwaysInline() void semidefiniteFactor(View a, Size size, int* order, MatrixView<T> factor)
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
      detail::exchange(a, n, k, pivot);
      const int swapped = order[k];
      order[k] = order[pivot];
      order[pivot] = swapped;
    }
    const T d = a(k, k);
    for (int i = k + 1; i < n; ++i) {
      a(i, k) /= d;
    }
    for (int i = k + 1; i < n; ++i) {
      for (int j = k + 1; j <= i; ++j) {
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

#endif // GTRK_BASE_DENSESOLVE_H_
