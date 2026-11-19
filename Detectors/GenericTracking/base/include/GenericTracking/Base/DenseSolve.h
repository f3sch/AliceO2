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
///        positive semi-definite matrix; the packed index of a symmetric
///        matrix (`covIndex`).
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

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/ElementaryMath.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

#include <type_traits>

namespace o2::gtrk
{

/// Index of element (row, col) of a symmetric matrix in its packed lower
/// triangle, row-major:
///
///     0
///     1  2
///     3  4  5
///     6  7  8  9
///    10 11 12 13 14
///
/// For the 5x5 track covariance this is the layout O2 uses: o2::track::CovarMap
/// is {{0,1,3,6,10},...} with DiagMap {0,2,5,9,14}, which is what this formula
/// produces.  O2 tabulates it and we compute it.
///
/// Symmetric in its arguments, so a caller never has to know which of the two
/// indices is the larger one -- which is exactly the mistake that makes a
/// packed covariance silently wrong rather than loudly broken.
GPUhdi() GPUAlwaysInline() constexpr int covIndex(int row, int col)
{
  const int hi = row > col ? row : col;
  const int lo = row > col ? col : row;
  return (hi * (hi + 1)) / 2 + lo;
}

/// The entries a packed symmetric n x n matrix holds, n (n + 1) / 2.
GPUhdi() GPUAlwaysInline() constexpr int covSize(int n) { return (n * (n + 1)) / 2; }

} // namespace o2::gtrk

namespace o2::gtrk::dense
{

/// A row-major matrix: element (row, column) at `data[row * stride + column]`.
template <typename T>
struct MatrixView {
  T* data = nullptr;
  int stride = 0;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() T& operator()(int row, int column) const { return data[(row * stride) + column]; }
};

/// A symmetric matrix by its lower triangle, packed (`covIndex`):
/// (row, column) and (column, row) are the same element.
template <typename T>
struct PackedView {
  T* data = nullptr;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() T& operator()(int row, int column) const { return data[covIndex(row, column)]; }
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
  for (int j = 0; j < k; ++j) {
    o2::gpu::CAMath::Swap(a(k, j), a(p, j));
  }
  o2::gpu::CAMath::Swap(a(k, k), a(p, p));
  for (int i = k + 1; i < p; ++i) {
    o2::gpu::CAMath::Swap(a(i, k), a(p, i));
  }
  for (int i = p + 1; i < n; ++i) {
    o2::gpu::CAMath::Swap(a(i, k), a(i, p));
  }
}

/// A 2 x 2 pivot block D = [[d11, d21], [d21, d22]] scaled by its off-diagonal
/// d21, which Bunch-Kaufman's choice makes the block's largest entry:
/// D = d21 [[e11, 1], [1, e22]], D^-1 = scale [[e22, -1], [-1, e11]] with
/// scale = 1 / (d21 (e11 e22 - 1)). No product of two entries is formed, so a
/// block of any magnitude neither overflows nor underflows.
template <typename T>
struct ScaledBlock {
  T e11;
  T e22;
  T reduced; ///< e11 e22 - 1: det D / d21^2
  T scale;
  GPUhdi() GPUAlwaysInline() ScaledBlock(T d11, T d21, T d22) : e11(d11 / d21), e22(d22 / d21), reduced((e11 * e22) - T(1)), scale((T(1) / d21) / reduced) {}
  /// (x1, x2) D^-1, in place.
  template <typename U>
  GPUhdi() GPUAlwaysInline() void solve(U& x1, U& x2) const
  {
    const U first = scale * ((e22 * x1) - x2);
    const U second = scale * ((e11 * x2) - x1);
    x1 = first;
    x2 = second;
  }
};
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
    // A non-finite entry would defeat every comparison below (and select a
    // 2 x 2 block past the last row): refused like a singular block.
    if (!isFiniteOf(diagonal) || !isFiniteOf(columnLargest) || !(diagonal != T(0) || columnLargest != T(0))) {
      return false;
    }
    int step = 1;
    int pivot = k;
    if (!(diagonal >= alpha * columnLargest)) {
      // columnLargest > 0 here, so largestRow > k: a 2 x 2 block has its row.
      T rowLargest = T(0);
      for (int j = k; j < n; ++j) {
        if (j != largestRow) {
          rowLargest = o2::gpu::GPUCommonMath::Max(rowLargest, o2::gpu::GPUCommonMath::Abs(j < largestRow ? a(largestRow, j) : a(j, largestRow)));
        }
      }
      // rowLargest >= columnLargest (row largestRow holds a(largestRow, k)),
      // so the ratio is at most one: no product of two entries.
      if (diagonal >= alpha * columnLargest * (columnLargest / rowLargest)) {
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
      const detail::ScaledBlock<T> block(a(k, k), a(k + 1, k), a(k + 1, k + 1));
      if (!(block.reduced != T(0)) || !isFiniteOf(block.scale)) {
        return false;
      }
      for (int i = k + 2; i < n; ++i) {
        // (l_i1, l_i2) = (a_ik, a_ik+1) D^-1
        T first = a(i, k);
        T second = a(i, k + 1);
        block.solve(first, second);
        for (int j = k + 2; j <= i; ++j) {
          a(i, j) -= (first * a(j, k)) + (second * a(j, k + 1));
        }
      }
      for (int i = k + 2; i < n; ++i) {
        T first = a(i, k);
        T second = a(i, k + 1);
        block.solve(first, second);
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
      o2::gpu::CAMath::Swap(x[k], x[swaps[k]]);
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
      const detail::ScaledBlock<ElementOf<View>> block(ldl(k, k), ldl(k + 1, k), ldl(k + 1, k + 1));
      T first = x[k];
      T second = x[k + 1];
      block.solve(first, second);
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
      o2::gpu::CAMath::Swap(x[k], x[swaps[k]]);
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
      // det D = d21^2 (e11 e22 - 1), its logarithm without the product
      const T d21 = o2::gpu::GPUCommonMath::Abs(ldl(k + 1, k));
      const detail::ScaledBlock<T> block(ldl(k, k), ldl(k + 1, k), ldl(k + 1, k + 1));
      sum += (T(2) * logOf(d21)) + logOf(o2::gpu::GPUCommonMath::Abs(block.reduced));
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
      o2::gpu::CAMath::Swap(order[k], order[pivot]);
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
