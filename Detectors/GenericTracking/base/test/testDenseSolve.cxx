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

/// \file testDenseSolve.cxx
/// \brief The dense factorisations against matrices whose answers are known.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK DenseSolve
#include "GenericTracking/Base/DenseSolve.h"

#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk::dense
{
namespace
{

using Matrix = std::vector<double>;

MatrixView<double> view(Matrix& a, std::int32_t n) { return {.data = a.data(), .stride = n}; }

// The index of row i, column j of a row-major matrix of n columns.
std::uint32_t entry(std::int32_t n, std::int32_t i, std::int32_t j) { return static_cast<std::uint32_t>((i * n) + j); }

// The index of row i, column j <= i of a packed lower triangle (covIndex).
std::uint32_t packedEntry(std::int32_t i, std::int32_t j) { return static_cast<std::uint32_t>(covIndex(i, j)); }

// a b^T + c I for the rows of b (n x k): positive definite for c > 0,
// semi-definite of rank k for c = 0.
Matrix gram(const std::vector<std::vector<double>>& b, double c)
{
  const std::size_t n = b.size();
  Matrix a(n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      double sum = i == j ? c : 0.;
      for (std::size_t k = 0; k < b[i].size(); ++k) {
        sum += b[i][k] * b[j][k];
      }
      a[(i * n) + j] = sum;
    }
  }
  return a;
}

Matrix product(const Matrix& a, const Matrix& b, std::int32_t n)
{
  Matrix c(static_cast<std::size_t>(n * n), 0.);
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j < n; ++j) {
      for (std::int32_t k = 0; k < n; ++k) {
        c[entry(n, i, j)] += a[entry(n, i, k)] * b[entry(n, k, j)];
      }
    }
  }
  return c;
}

void requireIdentity(const Matrix& a, std::int32_t n, double tolerance)
{
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j < n; ++j) {
      BOOST_TEST(std::abs(a[entry(n, i, j)] - (i == j ? 1. : 0.)) < tolerance);
    }
  }
}

const std::vector<std::vector<double>> kRows{{1., 2., 0.5}, {-0.5, 1., 3.}, {2., -1., 1.}, {0.25, 0.5, -2.}};

// The packed lower triangle of the symmetric \a a.
Matrix packed(const Matrix& a, std::int32_t n)
{
  Matrix p(static_cast<std::size_t>(covSize(n)));
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j <= i; ++j) {
      p[packedEntry(i, j)] = a[entry(n, i, j)];
    }
  }
  return p;
}

// A reference: Gaussian elimination with partial pivoting, solving a x = b
// and returning ln |det a|.
double referenceSolve(Matrix a, std::int32_t n, std::vector<double>& b)
{
  double logDet = 0.;
  for (std::int32_t k = 0; k < n; ++k) {
    std::int32_t p = k;
    for (std::int32_t i = k + 1; i < n; ++i) {
      if (std::abs(a[entry(n, i, k)]) > std::abs(a[entry(n, p, k)])) {
        p = i;
      }
    }
    for (std::int32_t j = 0; j < n; ++j) {
      std::swap(a[entry(n, k, j)], a[entry(n, p, j)]);
    }
    std::swap(b[static_cast<std::uint32_t>(k)], b[static_cast<std::uint32_t>(p)]);
    logDet += std::log(std::abs(a[entry(n, k, k)]));
    for (std::int32_t i = k + 1; i < n; ++i) {
      const double f = a[entry(n, i, k)] / a[entry(n, k, k)];
      for (std::int32_t j = k; j < n; ++j) {
        a[entry(n, i, j)] -= f * a[entry(n, k, j)];
      }
      b[static_cast<std::uint32_t>(i)] -= f * b[static_cast<std::uint32_t>(k)];
    }
  }
  for (std::int32_t k = n - 1; k >= 0; --k) {
    for (std::int32_t j = k + 1; j < n; ++j) {
      b[static_cast<std::uint32_t>(k)] -= a[entry(n, k, j)] * b[static_cast<std::uint32_t>(j)];
    }
    b[static_cast<std::uint32_t>(k)] /= a[entry(n, k, k)];
  }
  return logDet;
}

// Bunch-Kaufman on the packed \a a against the reference: the solution
// and ln |det|; true if any 2 x 2 block was taken.
bool checkLDLT(const Matrix& a, std::int32_t n, double tolerance)
{
  Matrix p = packed(a, n);
  std::vector<std::int32_t> swaps(static_cast<std::size_t>(n));
  std::vector<std::int32_t> blocks(static_cast<std::size_t>(n));
  BOOST_REQUIRE(ldltFactor(PackedView<double>{p.data()}, n, swaps.data(), blocks.data()));
  std::vector<double> x(static_cast<std::size_t>(n));
  for (std::int32_t i = 0; i < n; ++i) {
    x[static_cast<std::uint32_t>(i)] = 1. + (0.37 * i) - (0.05 * i * i);
  }
  std::vector<double> reference = x;
  const double logDet = referenceSolve(a, n, reference);
  ldltSolve(PackedView<double>{p.data()}, n, swaps.data(), blocks.data(), x.data());
  for (std::int32_t i = 0; i < n; ++i) {
    BOOST_TEST(std::abs(x[static_cast<std::uint32_t>(i)] - reference[static_cast<std::uint32_t>(i)]) < tolerance * (1. + std::abs(reference[static_cast<std::uint32_t>(i)])));
  }
  BOOST_TEST(std::abs(ldltLogAbsDeterminant(PackedView<double>{p.data()}, n, blocks.data()) - logDet) < tolerance * (1. + std::abs(logDet)));
  bool twoByTwo = false;
  for (const std::int32_t b : blocks) {
    twoByTwo = twoByTwo || b == 2;
  }
  return twoByTwo;
}

} // namespace

BOOST_AUTO_TEST_CASE(CholeskyInvertsSolvesAndGivesTheDeterminant)
{
  constexpr std::int32_t n = 4;
  const Matrix a = gram(kRows, 0.7);
  Matrix l = a;
  BOOST_REQUIRE(choleskyFactor(view(l, n), n));
  Matrix inverse(static_cast<std::size_t>(n * n));
  choleskyInverse(view(l, n), n, view(inverse, n));
  requireIdentity(product(a, inverse, n), n, 1e-12);

  // the determinant: the reference elimination's pivots
  std::vector<double> unused(n, 0.);
  BOOST_TEST(std::abs(choleskyLogDeterminant(view(l, n), n) - referenceSolve(a, n, unused)) < 1e-12);

  // the forward solve whitens: |L^-1 x|^2 = x^T a^-1 x
  std::array<double, n> x{0.3, -1.2, 2., 0.5};
  double expected = 0.;
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j < n; ++j) {
      expected += x[static_cast<std::uint32_t>(i)] * inverse[entry(n, i, j)] * x[static_cast<std::uint32_t>(j)];
    }
  }
  forwardSolve(view(l, n), n, x.data());
  double whitened = 0.;
  for (const double v : x) {
    whitened += v * v;
  }
  BOOST_TEST(std::abs(whitened - expected) < 1e-12 * expected);
}

BOOST_AUTO_TEST_CASE(CholeskyRefusesWhatIsNotPositiveDefinite)
{
  constexpr std::int32_t n = 4;
  Matrix semidefinite = gram(kRows, 0.);
  BOOST_TEST(!choleskyFactor(view(semidefinite, n), n));
  Matrix indefinite = gram(kRows, 0.7);
  indefinite[(2 * n) + 2] = -1.;
  BOOST_TEST(!choleskyFactor(view(indefinite, n), n));
}

BOOST_AUTO_TEST_CASE(OneCodeForAFixedAndARuntimeSize)
{
  constexpr std::int32_t n = 4;
  const Matrix a = gram(kRows, 0.7);
  Matrix runtime = a;
  Matrix fixed = a;
  BOOST_REQUIRE(choleskyFactor(view(runtime, n), n));
  BOOST_REQUIRE(choleskyFactor(view(fixed, n), std::integral_constant<std::int32_t, n>{}));
  BOOST_TEST(runtime == fixed);
}

BOOST_AUTO_TEST_CASE(BunchKaufmanSolvesSymmetricIndefiniteSystems)
{
  // The Lagrange system [[N, C^T], [C, 0]], N singular along x: the zero
  // diagonal needs a 2 x 2 pivot, and C fixes N's null direction.
  {
    constexpr std::int32_t n = 3;
    const Matrix system{0., 0., 1.,
                        0., 2., 1.,
                        1., 1., 0.};
    BOOST_TEST(checkLDLT(system, n, 1e-13));
  }
  // KKT-shaped: N = B B^T of rank 3 of 6, two rows of C fixing its null
  // directions, two more; and random symmetric indefinite matrices.
  {
    const std::vector<std::vector<double>> b{{1., 0.5, -1.}, {0.2, -1., 0.3}, {2., 1., 0.}, {0., 0.7, 1.5}, {-1., 0., 0.4}, {0.3, 0.3, 0.3}};
    const Matrix nMatrix = gram(b, 0.);
    const std::vector<std::vector<double>> c{{1., 0., 2., 0., -1., 0.5}, {0., 1., -1., 1., 0., 0.}, {0.5, 0.5, 0., 0., 1., 1.}, {1., -1., 0., 2., 0., 0.}};
    const std::int32_t n = 6;
    const auto m = static_cast<std::int32_t>(c.size());
    const std::int32_t size = n + m;
    Matrix system(static_cast<std::size_t>(size * size), 0.);
    for (std::int32_t i = 0; i < n; ++i) {
      for (std::int32_t j = 0; j < n; ++j) {
        system[entry(size, i, j)] = nMatrix[entry(n, i, j)];
      }
    }
    for (std::int32_t r = 0; r < m; ++r) {
      for (std::int32_t j = 0; j < n; ++j) {
        system[entry(size, n + r, j)] = c[static_cast<std::uint32_t>(r)][static_cast<std::uint32_t>(j)];
        system[entry(size, j, n + r)] = c[static_cast<std::uint32_t>(r)][static_cast<std::uint32_t>(j)];
      }
    }
    checkLDLT(system, size, 1e-11);
  }
  std::uint32_t state = 12345U;
  const auto next = [&]() {
    state = (1103515245U * state) + 12345U;
    return (static_cast<double>((state >> 8) & 0xffffU) / 32768.) - 1.;
  };
  bool anyTwoByTwo = false;
  for (std::int32_t trial = 0; trial < 50; ++trial) {
    const std::int32_t n = 2 + (trial % 9);
    Matrix a(static_cast<std::size_t>(n * n));
    for (std::int32_t i = 0; i < n; ++i) {
      for (std::int32_t j = 0; j <= i; ++j) {
        const double v = next();
        // small diagonals push it to 2 x 2 pivots
        a[entry(n, i, j)] = i == j ? 0.01 * v : v;
        a[entry(n, j, i)] = a[entry(n, i, j)];
      }
    }
    anyTwoByTwo = checkLDLT(a, n, 1e-9) || anyTwoByTwo;
  }
  BOOST_TEST(anyTwoByTwo);
}

BOOST_AUTO_TEST_CASE(BunchKaufmanRefusesASingularSystem)
{
  constexpr std::int32_t n = 3;
  // rank 2: the third row the sum of the first two
  Matrix singular = packed({1., 2., 3.,
                            2., 1., 3.,
                            3., 3., 6.},
                           n);
  std::array<std::int32_t, n> swaps{};
  std::array<std::int32_t, n> blocks{};
  BOOST_TEST(!ldltFactor(PackedView<double>{singular.data()}, n, swaps.data(), blocks.data()));
  Matrix zero(static_cast<std::size_t>(covSize(n)), 0.);
  BOOST_TEST(!ldltFactor(PackedView<double>{zero.data()}, n, swaps.data(), blocks.data()));
}

BOOST_AUTO_TEST_CASE(BunchKaufmanRefusesANonFiniteSystem)
{
  // A NaN pivot defeats every comparison: refused, nothing written beyond the
  // matrix's own rows (the sentinels after them stay).
  std::array<double, 1> nan{std::numeric_limits<double>::quiet_NaN()};
  std::array<std::int32_t, 2> swaps{-7, -7};
  std::array<std::int32_t, 2> blocks{-7, -7};
  BOOST_TEST(!ldltFactor(PackedView<double>{nan.data()}, 1, swaps.data(), blocks.data()));
  BOOST_TEST(swaps[1] == -7);
  BOOST_TEST(blocks[1] == -7);
  Vec3D infinite{1., std::numeric_limits<double>::infinity(), 2.};
  BOOST_TEST(!ldltFactor(PackedView<double>{infinite.data()}, 2, swaps.data(), blocks.data()));
}

BOOST_AUTO_TEST_CASE(BunchKaufmanHoldsAnyScale)
{
  // [[0, t], [t, 0]] x = [1, 1] is x = [1/t, 1/t] and ln |det| = 2 ln t, for
  // any t the type represents with its inverse: the 2 x 2 block's arithmetic
  // forms no product of two entries.
  for (const double t : {1e-200, 1e-100, 1., 1e100, 1e200}) {
    Vec3D a{0., t, 0.};
    std::array<std::int32_t, 2> swaps{};
    std::array<std::int32_t, 2> blocks{};
    BOOST_REQUIRE(ldltFactor(PackedView<double>{a.data()}, 2, swaps.data(), blocks.data()));
    Vec2D x{1., 1.};
    ldltSolve(PackedView<double>{a.data()}, 2, swaps.data(), blocks.data(), x.data());
    BOOST_TEST(std::abs((x[0] * t) - 1.) < 1e-15);
    BOOST_TEST(std::abs((x[1] * t) - 1.) < 1e-15);
    BOOST_TEST(std::abs(ldltLogAbsDeterminant(PackedView<double>{a.data()}, 2, blocks.data()) - (2. * std::log(t))) < 1e-12 * (1. + std::abs(2. * std::log(t))));
  }
}

BOOST_AUTO_TEST_CASE(PackedAndRowMajorCholeskyAgree)
{
  constexpr std::int32_t n = 4;
  const Matrix a = gram(kRows, 0.7);
  Matrix full = a;
  Matrix p = packed(a, n);
  BOOST_REQUIRE(choleskyFactor(view(full, n), n));
  BOOST_REQUIRE(choleskyFactor(PackedView<double>{p.data()}, n));
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j <= i; ++j) {
      BOOST_TEST(full[entry(n, i, j)] == p[packedEntry(i, j)]);
    }
  }
  Matrix inverse(static_cast<std::size_t>(covSize(n)));
  choleskyInverse(PackedView<double>{p.data()}, n, PackedView<double>{inverse.data()});
  Matrix unpacked(static_cast<std::size_t>(n * n));
  for (std::int32_t i = 0; i < n; ++i) {
    for (std::int32_t j = 0; j < n; ++j) {
      unpacked[entry(n, i, j)] = inverse[packedEntry(i, j)];
    }
  }
  requireIdentity(product(a, unpacked, n), n, 1e-12);
}

BOOST_AUTO_TEST_CASE(SemidefiniteFactorRebuildsTheMatrix)
{
  // rank 2 of 4 (its elimination exact in binary, so what is left is exactly
  // zero), and rank 4 (definite): F F^T = a, the columns past the rank zero
  for (const double c : {0., 0.5}) {
    constexpr std::int32_t n = 4;
    const Matrix a = gram({{1., 0.}, {0., 1.}, {1., 1.}, {1., -1.}}, c);
    Matrix work = packed(a, n);
    std::array<std::int32_t, n> order{};
    Matrix factor(static_cast<std::size_t>(n * n));
    semidefiniteFactor(PackedView<double>{work.data()}, n, order.data(), view(factor, n));
    for (std::int32_t i = 0; i < n; ++i) {
      for (std::int32_t j = 0; j < n; ++j) {
        double sum = 0.;
        for (std::int32_t k = 0; k < n; ++k) {
          sum += factor[entry(n, i, k)] * factor[entry(n, j, k)];
        }
        BOOST_TEST(std::abs(sum - a[entry(n, i, j)]) < 1e-12);
      }
    }
    if (exactlyEqual(c, 0.)) {
      for (std::int32_t i = 0; i < n; ++i) {
        for (std::int32_t k = 2; k < n; ++k) {
          BOOST_TEST(factor[entry(n, i, k)] == 0.);
        }
      }
    }
  }
  // a zero matrix: a zero factor
  constexpr std::int32_t n = 3;
  Matrix zero(static_cast<std::size_t>(n * n), 0.);
  std::array<std::int32_t, n> order{};
  Matrix factor(static_cast<std::size_t>(n * n), 1.);
  semidefiniteFactor(view(zero, n), n, order.data(), view(factor, n));
  for (const double v : factor) {
    BOOST_TEST(v == 0.);
  }
}

} // namespace o2::gtrk::dense
