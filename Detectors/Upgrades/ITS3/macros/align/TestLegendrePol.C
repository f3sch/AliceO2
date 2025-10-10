// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "TRandom.h"

#include "MathUtils/LegendrePols.h"
#endif

void TestLegendrePol()
{
  constexpr int nMaxOrder{2};
  constexpr int nPoints{100};
  constexpr int nPoints2{nPoints * nPoints};
  constexpr double minX{-1.0}, maxX{1.0},
    stepX{(maxX - minX) / static_cast<double>(nPoints)};

  gRandom->SetSeed(0);
  auto getRandom = []() {
    constexpr double scale{80.e-4};
    return scale * gRandom->Uniform(-1.0, 1.0);
  };

  { // 1D
    Info("", "---------------- 1D -------------");
    std::array<double, nPoints> x, y;
    for (int nOrder{0}; nOrder <= nMaxOrder; ++nOrder) {
      std::vector<double> coeff(nOrder + 1, 0.0);
      std::generate(std::begin(coeff), std::end(coeff), getRandom);

      o2::math_utils::Legendre1DPolynominal leg1D(coeff);

      auto c1d = new TCanvas(Form("c1D_%d", nOrder),
                             Form("Legendre 1D Order %d", nOrder));
      c1d->Divide(2, 1);

      { // Draw total polynominal
        for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
          double xcur = minX + iPoint * stepX;
          x[iPoint] = xcur;
          y[iPoint] = leg1D(xcur);
        }
        c1d->cd(2);
        auto g = new TGraph(nPoints, x.data(), y.data());
        g->SetName(Form("g1d_%d", nOrder));
        g->SetTitle(Form("Legendre Polynominal %d-deg;u;w", nOrder));
        g->Draw();
      }

      { // Draw single contributions
        auto mg = new TMultiGraph();
        auto leg = new TLegend();
        for (int iOrder{0}; iOrder <= nOrder; ++iOrder) {
          for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
            double xcur = minX + iPoint * stepX;
            x[iPoint] = xcur;
            y[iPoint] = leg1D(iOrder, xcur);
          }
          auto g = new TGraph(nPoints, x.data(), y.data());
          g->SetName(Form("g1d_%d_%d", nOrder, iOrder));
          g->SetTitle(Form("c_{%d}=%.4f;u;w", iOrder, coeff[iOrder]));
          mg->Add(g, "PL");
          leg->AddEntry(g, Form("c_{%d}=%.4f", iOrder, coeff[iOrder]), "lp");
        }
        c1d->cd(1);
        mg->SetTitle("Contribution decomposition;u;w");
        mg->Draw("A pmc plc");
        leg->Draw();
        gPad->Update();
      }
    }
  }

  { // 2D
    Info("", "---------------- 2D -------------");
    std::array<double, nPoints2> x, y, z;
    for (int nOrder{0}; nOrder <= nMaxOrder; ++nOrder) {
      auto c2d = new TCanvas(Form("c2D_%d", nOrder),
                             Form("Legendre 2D Order %d", nOrder));
      c2d->Divide(nOrder + 1, nOrder + 1);

      TMatrixD coeff(nOrder + 1, nOrder + 1);
      // Random initialization
      for (int i{0}; i <= nOrder; ++i) {
        for (int j{0}; j <= i; ++j) {
          coeff(i, j) = getRandom();
        }
      }

      o2::math_utils::Legendre2DPolynominal leg2D(coeff);
      leg2D.printCoefficients();

      { // Draw total polynominal
        for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
          double xcur = minX + iPoint * stepX;
          for (int jPoint{0}; jPoint < nPoints; ++jPoint) {
            double ycur = minX + jPoint * stepX;
            int i = iPoint * nPoints + jPoint;
            x[i] = xcur;
            y[i] = ycur;
            z[i] = leg2D(xcur, ycur);
          }
        }
        c2d->cd(nOrder + 1);
        auto g = new TGraph2D(nPoints2, x.data(), y.data(), z.data());
        g->SetTitle(Form("Legendre Polynominal %d-deg;u;v;w", nOrder));
        g->SetName(Form("g2d_%d", nOrder));
        g->Draw("surf1");
        g->GetXaxis()->SetRangeUser(minX, maxX);
        g->GetYaxis()->SetRangeUser(minX, maxX);
        gPad->Update();
      }

      { // Draw decomposition of contributions to polynominal
        for (int iOrder{0}; iOrder <= nOrder; ++iOrder) {
          for (int jOrder{0}; jOrder <= iOrder; ++jOrder) {
            for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
              double xcur = minX + iPoint * stepX;
              for (int jPoint{0}; jPoint < nPoints; ++jPoint) {
                double ycur = minX + jPoint * stepX;
                int i = iPoint * nPoints + jPoint;
                x[i] = xcur;
                y[i] = ycur;
                z[i] = leg2D(iOrder, jOrder, xcur, ycur);
              }
            }
            c2d->cd(1 + iOrder * (nOrder + 1) + jOrder);
            auto g = new TGraph2D(nPoints2, x.data(), y.data(), z.data());
            g->SetTitle(Form("c_{%d,%d}=%.4f;u;v;w", iOrder, jOrder,
                             coeff(iOrder, jOrder)));
            g->SetName(Form("g2d_%d_%d_%d", nOrder, iOrder, jOrder));
            if (iOrder == 0 && jOrder == 0) { // Fix display of constant value
              g->Draw("P0");
            } else {
              g->Draw("surf1");
            }
            g->GetXaxis()->SetRangeUser(minX, maxX);
            g->GetYaxis()->SetRangeUser(minX, maxX);
            gPad->Update();
          }
        }
      }
      c2d->Draw();
    }
  }
}
