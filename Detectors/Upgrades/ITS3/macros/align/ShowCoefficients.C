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
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "TRandom.h"
#include "TMath.h"

#include "MathUtils/LegendrePols.h"
#include "ITS3Align/Deformations.h"
#include "MathUtils/Utils.h"
#endif

static ROOT::Math::Minimizer* gMin;

void ShowCoefficients(const std::string& fileName = "misparams.root", bool findMin = false)
{
  o2::its3::align::Deformations def;
  def.init(fileName);

  if (findMin) {
    gMin = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if (gMin == nullptr) {
      Error("", "Cannot create minimizer !");
      return;
    }
    gMin->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    gMin->SetTolerance(0.00001);
    gMin->SetPrintLevel(1);
  }

  if (0) {
    const std::array<const char*, 3> axisName{"x", "y", "z"};
    constexpr int nPoints{100};
    constexpr int nPoints2{nPoints * nPoints};
    constexpr double minX{-1.0}, maxX{1.0},
      stepX{(maxX - minX) / static_cast<double>(nPoints)};

    for (unsigned int iSensor{0}; iSensor < 6; ++iSensor) {
      const auto nOrders = def.getOrders(iSensor);
      for (unsigned int iAxis{0}; iAxis < 3; ++iAxis) {
        std::array<double, nPoints2> x, y, z;
        auto canv = new TCanvas(Form("c_sensor%d_d%s", iSensor, axisName[iAxis]), Form("Legendre Coefficients Sensor %d - Axis %s", iSensor, axisName[iAxis]));
        canv->Divide(nOrders[iAxis] + 1, nOrders[iAxis] + 1);

        { // Draw total polynominal
          for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
            double xcur = minX + iPoint * stepX;
            for (int jPoint{0}; jPoint < nPoints; ++jPoint) {
              double ycur = minX + jPoint * stepX;
              int i = iPoint * nPoints + jPoint;
              x[i] = xcur;
              y[i] = ycur;
              z[i] = def.getDeformation(iSensor, iAxis, xcur, ycur);
            }
          }
          canv->cd(nOrders[iAxis] + 1);
          auto g = new TGraph2D(nPoints2, x.data(), y.data(), z.data());
          g->SetTitle(Form("Legendre Polynominal %d-deg Sensor %d #Delta%s;u;v;w", nOrders[iAxis], iSensor, axisName[iAxis]));
          g->SetName(Form("g_%d_%s", iSensor, axisName[iAxis]));
          g->Draw("surf1");
          g->GetXaxis()->SetRangeUser(minX, maxX);
          g->GetYaxis()->SetRangeUser(minX, maxX);
          gPad->Update();
        }

        { // Draw decomposition of contributions to polynominal
          const auto& leg = def.getLegendre(iSensor, iAxis);
          const auto coeff = leg.getCoefficients();
          for (unsigned int iOrder{0}; iOrder <= nOrders[iAxis]; ++iOrder) {
            for (unsigned int jOrder{0}; jOrder <= iOrder; ++jOrder) {
              for (int iPoint{0}; iPoint < nPoints; ++iPoint) {
                double xcur = minX + iPoint * stepX;
                for (int jPoint{0}; jPoint < nPoints; ++jPoint) {
                  double ycur = minX + jPoint * stepX;
                  int i = iPoint * nPoints + jPoint;
                  x[i] = xcur;
                  y[i] = ycur;
                  z[i] = leg(iOrder, jOrder, xcur, ycur);
                }
              }
              canv->cd(1 + iOrder * (nOrders[iAxis] + 1) + jOrder);
              auto g = new TGraph2D(nPoints2, x.data(), y.data(), z.data());
              g->SetTitle(Form("c_{%d,%d}=%.4f;u;v;w", iOrder, jOrder, coeff(iOrder, jOrder)));
              g->SetName(Form("g_%d_%d_%d_%d", iSensor, iAxis, iOrder, jOrder));
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

        canv->Draw();
        canv->SaveAs(Form("its3_sensor%d_%s.pdf", iSensor, axisName[iAxis]));
      }
    }
  }

  {
    constexpr int nPoints{50};
    constexpr int nPoints2{nPoints * nPoints};
    constexpr double radL = o2::its3::constants::radii[2] + 0.3, zL = o2::its3::constants::segment::lengthSensitive / 2.0 + 2.0;

    // Plot the whole geometry
    std::array<TGraph2D*, 6> gIdeal;
    std::array<TGraph2D*, 6> gDeformed;
    constexpr double z1{-o2::its3::constants::segment::lengthSensitive / 2.0}, z2{o2::its3::constants::segment::lengthSensitive / 2.0}, zTot{z2 - z1}, zStep{zTot / (nPoints - 1)};
    auto canv = new TCanvas();
    canv->Divide(2, 1);
    for (unsigned int iSensor{0}; iSensor < 6; ++iSensor) {
      std::array<double, nPoints2> xi, yi, zi, xd, yd, zd;
      const double radius = o2::its3::constants::radii[iSensor / 2];
      const bool isTop = iSensor % 2 == 0;
      const double phi1 = o2::math_utils::to02Pi(((isTop) ? 0.f : 1.f) * TMath::Pi() + std::asin(o2::its3::constants::equatorialGap / 2.f / radius));
      const double phi2 = o2::math_utils::to02Pi(((isTop) ? 1.f : 2.f) * TMath::Pi() - std::asin(o2::its3::constants::equatorialGap / 2.f / radius));
      const double phiTot{phi2 - phi1}, phiStep{phiTot / (nPoints - 1)};
      for (int iZ{0}; iZ < nPoints; ++iZ) {
        double z = z1 + iZ * zStep;
        for (int iPhi{0}; iPhi < nPoints; ++iPhi) {
          int i = iZ * nPoints + iPhi;
          double phi = phi1 + iPhi * phiStep;

          xi[i] = radius * std::cos(phi);
          yi[i] = radius * std::sin(phi);
          zi[i] = z;

          const double u = ((phi - phi1) * 2.f) / phiTot - 1.f;
          const double v = ((z - z1)) / zTot - 1.f;
          const auto [dx, dy, dz] = def.getDeformation(iSensor, u, v);
          xd[i] = xi[i] + dx;
          yd[i] = yi[i] + dy;
          zd[i] = zi[i] + dz;
        }
      }

      canv->cd(1);
      auto ideal = new TGraph2D(Form("g_ideal_%d", iSensor), "Ideal Geometry", nPoints2, xi.data(), zi.data(), yi.data());
      ideal->SetMarkerStyle(kFullCircle);
      ideal->SetMarkerSize(1);
      ideal->SetMarkerColor(kBlue);
      ideal->SetLineColor(kBlue);
      if (iSensor == 0) {
        ideal->Draw("p0");
      } else {
        ideal->Draw("p0;same");
        if (iSensor == 5) {
          gPad->Update();
          ideal->GetXaxis()->SetTitle("X");
          ideal->GetYaxis()->SetTitle("Z");
          ideal->GetZaxis()->SetTitle("Y");
          ideal->GetXaxis()->SetRangeUser(-radL, radL);
          ideal->GetZaxis()->SetRangeUser(-radL, radL);
          ideal->GetYaxis()->SetRangeUser(-zL, zL);
        }
      }

      canv->cd(2);
      auto deformed = new TGraph2D(Form("g_def_%d", iSensor), "Deformed Geometry", nPoints2, xd.data(), zd.data(), yd.data());
      deformed->SetMarkerStyle(kFullCircle);
      deformed->SetMarkerSize(1);
      deformed->SetMarkerColor(kRed);
      deformed->SetLineColor(kRed);
      if (iSensor == 0) {
        deformed->Draw("p0");
      } else {
        deformed->Draw("p0;same");
        if (iSensor == 5) {
          gPad->Update();
          deformed->GetXaxis()->SetTitle("X");
          deformed->GetYaxis()->SetTitle("Z");
          deformed->GetZaxis()->SetTitle("Y");
          deformed->GetXaxis()->SetRangeUser(-radL, radL);
          deformed->GetZaxis()->SetRangeUser(-radL, radL);
          deformed->GetYaxis()->SetRangeUser(-zL, zL);
        }
      }
    }

    // Optionally find a deformation
    /* if (findMin2D) { */
    /*   std::vector<double> cccc(nOrder + 2, 0.0); */
    /*   cccc[0] = nOrder; */
    /*   for (int i{0}; i <= nOrder; ++i) { */
    /*     for (int j{0}; j <= i; ++j) { */
    /*       int k = i * (i + 1) / 2 + j; */
    /*       cccc[1 + k] = coeff(i, j); */
    /*     } */
    /*   } */
    /*   auto ig = legendre_poly2D_integral(cccc.data()); */

    /*   gMin->Clear(); */
    /*   ROOT::Math::Functor fmin(&legendre_poly2D_integral, */
    /*                            2 + nOrder * (nOrder + 1) / 2 + nOrder); */
    /*   Info("", "ig=%f    parameters=%d", ig, */
    /*        2 + nOrder * (nOrder + 1) / 2 + nOrder); */
    /*   gMin->SetFunction(fmin); */
    /*   constexpr double minStep{0.001}; */
    /*   gMin->SetFixedVariable(0, "nOrder", nOrder); */
    /*   gMin->SetFixedVariable(1, "c_00", coeff(0, 0)); */
    /*   for (int iOrder{1}; iOrder <= nOrder; ++iOrder) { */
    /*     for (int jOrder{0}; jOrder <= iOrder; ++jOrder) { */
    /*       int k = iOrder * (iOrder + 1) / 2 + jOrder + 1; */
    /*       Info("", "Setting parameter %d", k); */
    /*       if (getRandom() < 0.0) { */
    /*         gMin->SetFixedVariable(k, Form("c_%d_%d", iOrder, jOrder), */
    /*                                coeff(iOrder, jOrder)); */
    /*       } else { */
    /*         gMin->SetVariable(k, Form("c_%d_%d", iOrder, jOrder), */
    /*                           coeff(iOrder, jOrder), minStep); */
    /*       } */
    /*     } */
    /*   } */
    /*   gMin->Minimize(); */
    /*   return; */
    /*   auto stat = gMin->Status(); */
    /*   auto min = gMin->MinValue(); */
    /*   if ((stat == 0 || stat == 1) && min < 0.01) { */
    /*     Info("", "Minimizer converged with %f; using new values!", min); */
    /*     const double *cmins = gMin->X(); */
    /*     for (int iOrder{1}; iOrder <= nOrder; ++iOrder) { */
    /*       for (int jOrder{0}; jOrder <= iOrder; ++jOrder) { */
    /*         int k = iOrder * (iOrder + 1) / 2 + jOrder; */
    /*         coeff(iOrder, jOrder) = cmins[k + 1]; */
    /*       } */
    /*     } */
    /*   } */
    /* } */
  }
}
