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
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/DenseSolve.h"

#include <array>
#include <cmath>
#include <type_traits>
#include <vector>

namespace o2::gtrk::dense
{
namespace
{

using Matrix = std::vector<double>;

MatrixView<double> view(Matrix& a, int n) { return {a.data(), n}; }

// a b^T + c I for the rows of b (n x k): positive definite for c > 0,
// semi-definite of rank k for c = 0.
Matrix gram(const std::vector<std::vector<double>>& b, double c)
{
  const int n = static_cast<int>(b.size());
  Matrix a(static_cast<std::size_t>(n * n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = i == j ? c : 0.;
      for (std::size_t k = 0; k < b[i].size(); ++k) {
        sum += b[i][k] * b[j][k];
      }
      a[(i * n) + j] = sum;
    }
  }
  return a;
}

Matrix product(const Matrix& a, const Matrix& b, int n)
{
  Matrix c(static_cast<std::size_t>(n * n), 0.);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        c[(i * n) + j] += a[(i * n) + k] * b[(k * n) + j];
      }
    }
  }
  return c;
}

void requireIdentity(const Matrix& a, int n, double tolerance)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      BOOST_TEST(std::abs(a[(i * n) + j] - (i == j ? 1. : 0.)) < tolerance);
    }
  }
}

const std::vector<std::vector<double>> kRows{{1., 2., 0.5}, {-0.5, 1., 3.}, {2., -1., 1.}, {0.25, 0.5, -2.}};

} // namespace

BOOST_AUTO_TEST_CASE(CholeskyInvertsSolvesAndGivesTheDeterminant)
{
  constexpr int n = 4;
  const Matrix a = gram(kRows, 0.7);
  Matrix l = a;
  BOOST_REQUIRE(choleskyFactor(view(l, n), n));
  Matrix inverse(static_cast<std::size_t>(n * n));
  choleskyInverse(view(l, n), n, view(inverse, n));
  requireIdentity(product(a, inverse, n), n, 1e-12);

  // the determinant: the product of LU's pivots of the same matrix
  Matrix lu = a;
  std::array<int, n> pivots{};
  BOOST_REQUIRE(luFactor(view(lu, n), n, pivots.data()));
  BOOST_TEST(std::abs(choleskyLogDeterminant(view(l, n), n) - luLogAbsDeterminant(view(lu, n), n)) < 1e-12);

  // the forward solve whitens: |L^-1 x|^2 = x^T a^-1 x
  std::array<double, n> x{0.3, -1.2, 2., 0.5};
  double expected = 0.;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      expected += x[i] * inverse[(i * n) + j] * x[j];
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
  constexpr int n = 4;
  Matrix semidefinite = gram(kRows, 0.);
  BOOST_TEST(!choleskyFactor(view(semidefinite, n), n));
  Matrix indefinite = gram(kRows, 0.7);
  indefinite[(2 * n) + 2] = -1.;
  BOOST_TEST(!choleskyFactor(view(indefinite, n), n));
}

BOOST_AUTO_TEST_CASE(OneCodeForAFixedAndARuntimeSize)
{
  constexpr int n = 4;
  const Matrix a = gram(kRows, 0.7);
  Matrix runtime = a;
  Matrix fixed = a;
  BOOST_REQUIRE(choleskyFactor(view(runtime, n), n));
  BOOST_REQUIRE(choleskyFactor(view(fixed, n), std::integral_constant<int, n>{}));
  BOOST_TEST(runtime == fixed);
}

BOOST_AUTO_TEST_CASE(LUPivotsAndSolvesAnIndefiniteSystem)
{
  // The Lagrange system [[N, C^T], [C, 0]], N singular along x: the zero
  // diagonal needs row exchanges, and C fixes N's null direction.
  constexpr int n = 3;
  Matrix system{0., 0., 1.,
                0., 2., 1.,
                1., 1., 0.};
  const Matrix original = system;
  std::array<int, n> pivots{};
  BOOST_REQUIRE(luFactor(view(system, n), n, pivots.data()));
  BOOST_TEST(pivots[0] != 0);
  std::array<double, n> x{1., 4., 3.};
  const auto rhs = x;
  luSolve(view(system, n), n, pivots.data(), x.data());
  for (int i = 0; i < n; ++i) {
    double sum = 0.;
    for (int j = 0; j < n; ++j) {
      sum += original[(i * n) + j] * x[j];
    }
    BOOST_TEST(std::abs(sum - rhs[i]) < 1e-12);
  }
  // det = -2: |det| from the pivots
  BOOST_TEST(std::abs(luLogAbsDeterminant(view(system, n), n) - std::log(2.)) < 1e-12);
}

BOOST_AUTO_TEST_CASE(LURefusesASingularSystem)
{
  constexpr int n = 3;
  Matrix singular{1., 2., 3.,
                  2., 4., 6.,
                  0., 1., 1.};
  std::array<int, n> pivots{};
  BOOST_TEST(!luFactor(view(singular, n), n, pivots.data()));
}

BOOST_AUTO_TEST_CASE(SemidefiniteFactorRebuildsTheMatrix)
{
  // rank 2 of 4 (its elimination exact in binary, so what is left is exactly
  // zero), and rank 4 (definite): F F^T = a, the columns past the rank zero
  for (const double c : {0., 0.5}) {
    constexpr int n = 4;
    const Matrix a = gram({{1., 0.}, {0., 1.}, {1., 1.}, {1., -1.}}, c);
    Matrix work = a;
    std::array<int, n> order{};
    Matrix factor(static_cast<std::size_t>(n * n));
    semidefiniteFactor(view(work, n), n, order.data(), view(factor, n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double sum = 0.;
        for (int k = 0; k < n; ++k) {
          sum += factor[(i * n) + k] * factor[(j * n) + k];
        }
        BOOST_TEST(std::abs(sum - a[(i * n) + j]) < 1e-12);
      }
    }
    if (c == 0.) {
      for (int i = 0; i < n; ++i) {
        for (int k = 2; k < n; ++k) {
          BOOST_TEST(factor[(i * n) + k] == 0.);
        }
      }
    }
  }
  // a zero matrix: a zero factor
  constexpr int n = 3;
  Matrix zero(static_cast<std::size_t>(n * n), 0.);
  std::array<int, n> order{};
  Matrix factor(static_cast<std::size_t>(n * n), 1.);
  semidefiniteFactor(view(zero, n), n, order.data(), view(factor, n));
  for (const double v : factor) {
    BOOST_TEST(v == 0.);
  }
}

} // namespace o2::gtrk::dense
