#if !defined(__CLING__) || defined(__ROOTCLING__)
// #define ENABLE_UPGRADES
#include "CCDB/CcdbApi.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include <TFile.h>
#include <TRandom.h>
#include <fmt/format.h>
#include <vector>
#endif

using AlgPar = std::array<double, 6>;

AlgPar generateMisalignment(double x, double y, double z, double psi,
                            double theta, double phi);

void ITS3Misaligner(const std::string &ccdbHost = "http://localhost:8080",
                    long tmin = 0, long tmax = -1, double xEnv = 0.,
                    double yEnv = 0., double zEnv = 0., double psiEnv = 0.,
                    double thetaEnv = 0., double phiEnv = 0., double xHBa = 0.,
                    double yHBa = 0., double zHBa = 0., double psiHBa = 0.,
                    double thetaHBa = 0., double phiHBa = 0.,
                    const std::string &objectPath = "",
                    const std::string &fileName = "ITSAlignment.root") {
  std::vector<o2::detectors::AlignParam> params;
  o2::base::GeometryManager::loadGeometry("", false);
  auto geom = o2::its::GeometryTGeo::Instance();
  std::string symname;
  AlgPar pars;
  bool glo = true;
  symname = geom->composeSymNameITS();

  o2::detectors::DetID detITS("IT3");

  // ITS envelope
  pars = generateMisalignment(xEnv, yEnv, zEnv, psiEnv, thetaEnv, phiEnv);
  params.emplace_back(symname.c_str(), -1, pars[0], pars[1], pars[2], pars[3],
                      pars[4], pars[5], glo);
      LOGP(info, "Adding for {} -> {} {} {} {} {} {}", symname, pars[0],
           pars[1], pars[2], pars[3], pars[4], pars[5]);

  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {

    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      symname = geom->composeSymNameHalfBarrel(ilr, ihb, ilr < 3);
      if (ilr == 1) {
        pars = AlgPar{-50e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
      } else {
        pars = generateMisalignment(xHBa, yHBa, zHBa, psiHBa, thetaHBa, phiHBa);
      }
      params.emplace_back(symname.c_str(), -1, pars[0], pars[1], pars[2],
                          pars[3], pars[4], pars[5], glo);
      LOGP(info, "Adding for {} -> {} {} {} {} {} {}", symname, pars[0],
           pars[1], pars[2], pars[3], pars[4], pars[5]);
    }
  }

  if (!ccdbHost.empty()) {
    std::string path =
        objectPath.empty()
            ? o2::base::DetectorNameConf::getAlignmentPath(detITS)
            : objectPath;
    LOGP(info, "Storing alignment object on {}/{}", ccdbHost, path);
    o2::ccdb::CcdbApi api;
    map<string, string> metadata; // can be empty
    api.init(
        ccdbHost.c_str()); // or http://localhost:8080 for a local installation
    // store abitrary user object in strongly typed manner
    api.storeAsTFileAny(&params, path, metadata, tmin, tmax);
  }
  if (!fileName.empty()) {
    LOGP(info, "Storing ITS3 alignment in local file {}", fileName);
    TFile algFile(fileName.c_str(), "recreate");
    algFile.WriteObjectAny(&params, "std::vector<o2::detectors::AlignParam>", "alignment");
    algFile.Close();
  }
}

AlgPar generateMisalignment(double x, double y, double z, double psi,
                            double theta, double phi) {
  AlgPar pars;
  pars[0] = gRandom->Gaus(0, x);
  pars[1] = gRandom->Gaus(0, y);
  pars[2] = gRandom->Gaus(0, z);
  pars[3] = gRandom->Gaus(0, psi);
  pars[4] = gRandom->Gaus(0, theta);
  pars[5] = gRandom->Gaus(0, phi);
  return std::move(pars);
}
