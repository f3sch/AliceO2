// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file LegendrePols.h
/// \brief Definition of the NDim Legendre Polynominals
/// \author felix.schlepper@cern.ch

#ifndef LEGENDRE_NDIM_POLYNOMINAL_H_
#define LEGENDRE_NDIM_POLYNOMINAL_H_

#include "TNamed.h"
#include "Math/IParamFunction.h"
#include "TMatrixD.h"

#include <vector>

#include <boost/math/special_functions/legendre.hpp>

namespace o2::math_utils
{

// Defines the 1D Legendre Polynominals with coefficients:
// w(u) = c_0 + c_1 * u + c_2 * 0.5 * (3 * u * u - 1) + ...
// for u in [-1.0, 1.0]
class Legendre1DPolynominal : public TNamed,
                              public ROOT::Math::IParametricFunctionOneDim
{
 public:
  Legendre1DPolynominal(unsigned int order) : fOrder(order) {}
  Legendre1DPolynominal(const std::vector<double> p)
    : fOrder(p.size() - 1), fParams(p) {}

  double operator()(double x) const { return DoEvalPar(x, Parameters()); }
  double operator()(int i, double x) const
  {
    return DoEvalParSingle(i, x, Parameters());
  }

  virtual const double* Parameters() const override { return &fParams.front(); }

  virtual void SetParameters(const double* p) override
  {
    fParams = std::vector<double>(p, p + NPar());
  }

  virtual unsigned int NPar() const override { return fParams.size(); }

  // Unimplemented
  virtual ROOT::Math::IBaseFunctionOneDim* Clone() const override { return nullptr; }
  virtual TObject* Clone(const char*) const override { return nullptr; }

 private:
  virtual double DoEvalPar(double x, const double* p) const override
  {
    double sum{0.0};
    for (unsigned int iOrder{0}; iOrder <= fOrder; ++iOrder) {
      sum += p[iOrder] * boost::math::legendre_p(iOrder, x);
    }
    return sum;
  }

  double DoEvalParSingle(int i, double x, const double* p) const
  {
    return p[i] * boost::math::legendre_p(i, x);
  }

  unsigned int fOrder{0};
  std::vector<double> fParams;

  ClassDefOverride(o2::math_utils::Legendre1DPolynominal, 1);
};

// Defines the 2D Legendre Polynominals with coefficients:
// w(u, v) = c_00 +
//           c_10 * u + c_11 * v +
//           c_20 * 0.5 * (3 * u * u - 1) + c_21 * u * v + c_22 * (3 * v * v - 1) +
///          ....
// for u&v in [-1.0, 1.0]
class Legendre2DPolynominal : public TNamed,
                              public ROOT::Math::IParametricFunctionMultiDim
{
 public:
  Legendre2DPolynominal(unsigned int order) : fOrder(order) {}
  Legendre2DPolynominal(const std::vector<double>& p)
    : fOrder(p.size() - 1), fParams(p) {}
  Legendre2DPolynominal(const TMatrixD& p) : fOrder(p.GetNrows() - 1)
  {
    fParams = std::vector<double>(NPar());
    for (unsigned int iOrder{0}; iOrder <= fOrder; ++iOrder) {
      for (unsigned int jOrder{0}; jOrder <= iOrder; ++jOrder) {
        fParams[getFlatIdx(iOrder, jOrder)] = p(iOrder, jOrder);
      }
    }
  }

  double operator()(const double* x) const
  {
    return DoEvalPar(x, Parameters());
  }
  double operator()(double x, double y) const { return DoEvalPar(x, y); }
  double operator()(int i, int j, const double* x) const
  {
    return DoEvalParSingle(i, j, x, Parameters());
  }
  double operator()(int i, int j, double x, double y) const
  {
    return DoEvalParSingle(i, j, x, y, Parameters());
  }

  virtual const double* Parameters() const override { return &fParams.front(); }

  virtual void SetParameters(const double* p) override
  {
    fParams = std::vector<double>(p, p + NPar());
  }

  virtual unsigned int NPar() const override
  {
    return fOrder * (fOrder + 1) / 2 + fOrder + 1;
  }
  virtual unsigned int NDim() const override { return 2; }

  TMatrixD getCoefficients() const
  {
    TMatrixD m(fOrder + 1, fOrder + 1);
    for (unsigned int iOrder{0}; iOrder <= fOrder; ++iOrder) {
      for (unsigned int jOrder{0}; jOrder <= iOrder; ++jOrder) {
        m(iOrder, jOrder) = fParams[getFlatIdx(iOrder, jOrder)];
      }
    }
    return m;
  }

  void printCoefficients() const { getCoefficients().Print(); }

  // Unimplemented
  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return nullptr; }
  virtual TObject* Clone(const char*) const override { return nullptr; }

 private:
  virtual double DoEvalPar(const double* x, const double* p) const override
  {
    double sum{0.0};
    for (unsigned int iOrder{0}; iOrder <= fOrder; ++iOrder) {
      for (unsigned int jOrder{0}; jOrder <= iOrder; ++jOrder) {
        sum += DoEvalParSingle(iOrder, jOrder, x, p);
      }
    }
    return sum;
  }

  double DoEvalPar(double x, double y) const
  {
    double sum{0.0};
    for (unsigned int iOrder{0}; iOrder <= fOrder; ++iOrder) {
      for (unsigned int jOrder{0}; jOrder <= iOrder; ++jOrder) {
        sum += DoEvalParSingle(iOrder, jOrder, x, y, Parameters());
      }
    }
    return sum;
  }

  double DoEvalParSingle(int i, int j, const double* x, const double* p) const
  {
    return DoEvalParSingle(i, j, x[0], x[1], p);
  }

  double DoEvalParSingle(int i, int j, double x, double y,
                         const double* p) const
  {
    return p[getFlatIdx(i, j)] * boost::math::legendre_p(j, x) *
           boost::math::legendre_p(i - j, y);
  }

  inline int getFlatIdx(int i, int j) const { return i * (i - 1) / 2 + j; }

  unsigned int fOrder{0};
  std::vector<double> fParams;

  ClassDefOverride(o2::math_utils::Legendre2DPolynominal, 1);
};

} // namespace o2::math_utils

#endif
