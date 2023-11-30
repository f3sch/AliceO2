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

/// \file tile.C
/// \brief Prototype of ITS3Layer
/// \author felix.schlepper@cern.ch

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "Rtypes.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TGeoCompositeShape.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TMath.h"

#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TEvePointSet.h"
#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"

#include <iostream>
#include <fmt/format.h>

#include "CommonConstants/MathConstants.h"
#include "MathUtils/Cartesian.h"
#endif

namespace my
{
namespace its3
{
namespace material
{

enum Type : uint8_t {
  Vacuum = 0,
  Silicon,
  DeadZone,
  NMaterials,
};

static constexpr std::array<const char*, NMaterials> MaterialNames{
  "ITS3_VACUUM",
  "ITS3_SILICON",
  "ITS3_DEADZONE",
};

class Material
{
 public:
  // Getters
  const TGeoMaterial* getMaterial(Type type) const { return mMaterials[type]; }
  const TGeoMedium* getMedium(Type type) const { return mMediums[type]; }
  TGeoMaterial* getMaterial(Type type) { return mMaterials[type]; }
  TGeoMedium* getMedium(Type type) { return mMediums[type]; }

  // Print
  void print() const
  {
    for (int i{0}; i < NMaterials; ++i) {
      printf("Material: %s\n", MaterialNames[i]);
      mMaterials[i]->Print();
      mMediums[i]->Print();
    }
  }

  static Material& Instance()
  {
    static Material mat;
    return mat;
  }

 private:
  Material()
  {
    make(Vacuum, 0, 0, 0);
    make(Silicon, 26.98, 13, 2.7);
    make(DeadZone, 999, 12, 999);
  }
  Material(const Material&) = delete;
  Material& operator=(const Material&) = delete;
  ~Material() = default;

  void make(Type type, double aMass, double aNumber, double rho,
            double radlen = 0, double intlen = 0)
  {
    mMaterials[type] = new TGeoMaterial(MaterialNames[type], aMass, aNumber,
                                        rho, radlen, intlen);
    mMediums[type] =
      new TGeoMedium(MaterialNames[type], mCounter++, mMaterials[type]);
  }

  std::array<TGeoMaterial*, NMaterials> mMaterials{};
  std::array<TGeoMedium*, NMaterials> mMediums{};
  int mCounter{0};
};
} // namespace material

namespace constants
{
constexpr double PI = 3.14159274101257324e+00f;
constexpr double Rad2Deg = 180.f / PI;
constexpr double Deg2Rad = PI / 180.f;

constexpr double cm{1e+2}; // This is the default unit of TGeo so we use this as scale
constexpr double mu{1e-6 * cm};
constexpr double mm{1e-3 * cm};
namespace pixelarray
{
constexpr double length{9.197 * mm};
constexpr double width{3.571 * mm};
constexpr EColor color{kGreen};
constexpr unsigned int nCols{156};
constexpr unsigned int nRows{440};
constexpr unsigned int nPixels{nCols * nRows};
namespace pixel
{
constexpr double pitchCol{width / static_cast<double>(nCols)};
constexpr double pitchRow{length / static_cast<double>(nRows)};
} // namespace pixel
} // namespace pixelarray
namespace tile
{
namespace biasing
{
constexpr double length{0.06 * mm};
constexpr double width{3.571 * mm};
constexpr EColor color{kYellow};
static_assert(width == pixelarray::width);
} // namespace biasing
namespace powerswitches
{
constexpr double length{9.257 * mm};
constexpr double width{0.02 * mm};
constexpr double z{pixelarray::width};
constexpr EColor color{kBlue};
} // namespace powerswitches
namespace readout
{
constexpr double length{0.525 * mm};
constexpr double width{3.591 * mm};
constexpr EColor color{kMagenta};
static_assert(width == (biasing::width + powerswitches::width));
} // namespace readout
constexpr double width{readout::width};
constexpr double length{powerswitches::length + readout::length};
} // namespace tile
namespace rsu
{
namespace databackbone
{
constexpr double length{9.782 * mm};
constexpr double width{0.06 * mm};
constexpr EColor color{kRed};
} // namespace databackbone
constexpr double length{19.564 * mm};
constexpr double width{21.666 * mm};
} // namespace rsu
namespace segment
{
constexpr double length{rsu::length};
namespace lec
{
constexpr double length{segment::length};
constexpr double width{4.5 * mm};
constexpr EColor color{kCyan};
} // namespace lec
namespace rec
{
constexpr double length{segment::length};
constexpr double width{1.5 * mm};
constexpr EColor color{kCyan};
} // namespace rec
constexpr unsigned int nRSUs{12};
constexpr double width{nRSUs * rsu::width + lec::width + rec::width};
} // namespace segment
namespace carbonfoam
{
// TODO: Waiting for the further information from WP5(Corrado)
constexpr double longeronsHeight{2.0 * mm};                                // what is the height of the longerons?
constexpr double longeronsWidth{263 * mm};                                 // from blueprint
constexpr double HringWidth{6.0 * mm};                                     // from blueprint
constexpr double edgeBetwChipAndFoam{1.0 * mm};                            // from blueprint but not used cause forms are already overlapping
constexpr double gapBetwHringsLongerons{0.05 * mm};                        // from blueprint
constexpr std::array<int, 3> nHoles{11, 11, 11};                           // how many holes for each layer?
constexpr std::array<double, 3> radiusHoles{1.0 * mm, 1.0 * mm, 2.0 * mm}; // what is the radius of the holes for each layer?
constexpr EColor color{kGray};
} // namespace carbonfoam
constexpr unsigned int nLayers{3};
constexpr std::array<double, nLayers> radii{19 * mm, 25.2 * mm, 31.5 * mm}; // middle radius e.g. inner radius+thickness/2.
constexpr double equatorialGap{1 * mm};
constexpr std::array<unsigned int, nLayers> nSegments{3, 4, 5};
constexpr double thickness{50 * mu};
constexpr double effThickness{66 * mu};
} // namespace constants

class ITS3Layer
{
  // The hierarchy will be the following:
  // ITS2          ->       ITS3
  // ---------------------------------
  // Sensor                 PixelArray
  // Chip                   Tile
  // Module                 RSU
  // HalfStave              Segment
  // Stave                  Chip
  // HalfBarrel             CarbonForm
  // Layer                  Layer
 public:
  // Create one layer of ITS3 and attach it to the motherVolume.
  void createLayer(TGeoVolume* motherVolume, int layer = 0,
                   bool verbose = false)
  {
    if (layer > 2) {
      return;
    }
    mVerbose = verbose;
    mNLayer = layer;

    init();
    createPixelArray();
    createTile();
    createRSU();
    createSegment();
    createChip();
    createCarbonForm();
    createLayerImpl();

    // Add it to motherVolume
    // We have to shift it to the middle.
    auto zMove = new TGeoTranslation(0, 0, -constants::rsu::width * (constants::segment::nRSUs / 2.));
    motherVolume->AddNode(mLayer, 0, zMove);
  }

 private:
  void init()
  {
    // First we start by creating variables we are reusing a couple of times.
    mR = constants::radii[mNLayer];
    mRmin = mR - constants::thickness / 2.;
    mRmax = mR + constants::thickness / 2.;

    if (mVerbose) {
      std::cout << "INIT: Layer=" << mNLayer << "mR=" << mR << " Rmin=" << mRmin << " Rmax=" << mRmax << std::endl;
    }
  }

  void createPixelArray()
  {
    using namespace constants::pixelarray;
    // A pixel array is pure silicon and the sensitive part of our detector.
    // It will be segmented into a 440x144 matrix by the
    // SuperSegmentationAlpide.
    // Pixel Array is just a longer version of the biasing but starts in phi at
    // biasPhi2.
    double pixelArrayPhi1 =
      constants::tile::biasing::length / mRmin * constants::Rad2Deg;
    double pixelArrayPhi2 =
      length / mRmin * constants::Rad2Deg + pixelArrayPhi1;
    auto pixelArray = new TGeoTubeSeg(mRmin, mRmax, width / 2.,
                                      pixelArrayPhi1, pixelArrayPhi2);
    mPixelArray = new TGeoVolume(
      Form("pixelarray_%d", mNLayer) /* TODO change to correct name */,
      pixelArray,
      gGeoManager->GetMedium(material::MaterialNames[material::Silicon]));
    mPixelArray->SetLineColor(color);
    mPixelArray->RegisterYourself();
    if (mVerbose) {
      std::cout << "PixelArray:" << std::endl;
      mPixelArray->InspectShape();
      mPixelArray->InspectMaterial();
    }
  }

  void createTile()
  {
    using namespace constants::tile;
    // This functions creates a single Tile, which is the basic building block
    // of the chip. It consists of a pixelArray (sensitive area), biasing, power
    // switches and readout periphery (latter three are insensitive). We start
    // building the tile with the left upper edge of the biasing as center of
    // the tile’s z-coordinate axis.
    mTile = new TGeoVolumeAssembly(Form("tile_%d", mNLayer));
    mTile->VisibleDaughters();

    // Biasing
    auto zMoveBiasing = new TGeoTranslation(0, 0, +biasing::width / 2.);
    double biasPhi1 = 0;
    double biasPhi2 = biasing::length / mRmin * constants::Rad2Deg;
    auto biasing =
      new TGeoTubeSeg(mRmin, mRmax, biasing::width / 2, biasPhi1,
                      biasPhi2);
    auto biasingVol = new TGeoVolume(
      Form("biasing_%d", mNLayer), biasing,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    biasingVol->SetLineColor(biasing::color);
    biasingVol->RegisterYourself();
    if (mVerbose) {
      std::cout << "Biasing:" << std::endl;
      biasingVol->InspectShape();
      biasingVol->InspectMaterial();
    }
    mTile->AddNode(biasingVol, 0, zMoveBiasing);

    // Pixel Array is just a longer version of the biasing but starts in phi at
    // biasPhi2.
    mTile->AddNode(mPixelArray, 0, zMoveBiasing);

    // The readout periphery is also on top of the pixel array but extrudes on +z a bit e.g. is wider.
    auto zMoveReadout = new TGeoTranslation(0, 0, +readout::width / 2.);
    double readoutPhi1 =
      constants::pixelarray::length / mRmin * constants::Rad2Deg + biasPhi2;
    double readoutPhi2 =
      readout::length / mRmin * constants::Rad2Deg + readoutPhi1;
    auto readout = new TGeoTubeSeg(mRmin, mRmax, readout::width / 2,
                                   readoutPhi1, readoutPhi2);
    auto readoutVol = new TGeoVolume(
      Form("readout_%d", mNLayer), readout,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    readoutVol->SetLineColor(readout::color);
    readoutVol->RegisterYourself();
    if (mVerbose) {
      std::cout << "Readout:" << std::endl;
      readoutVol->InspectShape();
      readoutVol->InspectMaterial();
    }
    mTile->AddNode(readoutVol, 0, zMoveReadout);

    // Power Switches are on the side right side of the pixel array and biasing.
    auto zMovePowerSwitches = new TGeoTranslation(0, 0, +powerswitches::width / 2. + biasing::width);
    double powerPhi1 = 0;
    double powerPhi2 = powerswitches::length / mRmin * constants::Rad2Deg;
    auto powerSwitches = new TGeoTubeSeg(
      mRmin, mRmax, powerswitches::width / 2, powerPhi1, powerPhi2);
    auto powerSwitchesVol = new TGeoVolume(
      Form("powerswitches_%d", mNLayer), powerSwitches,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    powerSwitchesVol->SetLineColor(powerswitches::color);
    if (mVerbose) {
      std::cout << "PowerSwitches:" << std::endl;
      powerSwitchesVol->InspectShape();
      powerSwitchesVol->InspectMaterial();
    }
    mTile->AddNode(powerSwitchesVol, 0, zMovePowerSwitches);

    if (mSubstrate) {
      // Create the substrate layer at the back of the tile.
      // TODO
    }
  }

  void createRSU()
  {
    using namespace constants::rsu;
    // A Repeated Sensor Unit (RSU) is 12 Tiles + 4 Databackbones stichted together.
    mRSU = new TGeoVolumeAssembly(Form("rsu_%d", mNLayer));
    mRSU->VisibleDaughters();
    int nCopyRSU{0}, nCopyDB{0};

    // Create the DatabackBone
    // The Databackbone spans the whole phi of the tile.
    double dataBackbonePhi1 = 0;
    double dataBackbonePhi2 = databackbone::length / mRmin *
                              constants::Rad2Deg;
    auto dataBackbone = new TGeoTubeSeg(mRmin, mRmax, databackbone::width / 2.,
                                        dataBackbonePhi1,
                                        dataBackbonePhi2);
    auto dataBackboneVol = new TGeoVolume(
      Form("databackbone_%d", mNLayer), dataBackbone,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    dataBackboneVol->SetLineColor(databackbone::color);
    dataBackboneVol->RegisterYourself();
    if (mVerbose) {
      std::cout << "DataBackbone:" << std::endl;
      dataBackboneVol->InspectShape();
      dataBackboneVol->InspectMaterial();
    }

    // Lower Left
    auto zMoveLL1 = new TGeoTranslation(0, 0, constants::tile::width);
    auto zMoveLL2 = new TGeoTranslation(0, 0, constants::tile::width * 2.);
    auto zMoveLLDB = new TGeoTranslation(0, 0, -databackbone::width / 2.);
    // Lets attach the tiles to the QS.
    mRSU->AddNode(mTile, nCopyRSU++, nullptr);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveLL1);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveLL2);
    mRSU->AddNode(dataBackboneVol, nCopyDB++, zMoveLLDB);

    // Lower Right
    auto zMoveLR0 = new TGeoTranslation(0, 0, +width / 2.);
    auto zMoveLR1 = new TGeoTranslation(0, 0, constants::tile::width + width / 2.);
    auto zMoveLR2 = new TGeoTranslation(0, 0, constants::tile::width * 2. + width / 2.);
    auto zMoveLRDB = new TGeoTranslation(0, 0, -databackbone::width / 2. + width / 2.);
    // Lets attach the tiles to the QS.
    mRSU->AddNode(mTile, nCopyRSU++, zMoveLR0);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveLR1);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveLR2);
    mRSU->AddNode(dataBackboneVol, nCopyDB++, zMoveLRDB);

    // Rotation for top half
    double phi = length / mRmin * constants::Rad2Deg;
    auto rot = new TGeoRotation("", 0, 0, phi / 2.);

    // Upper Left
    auto zMoveUL1 = new TGeoCombiTrans(0, 0, constants::tile::width, rot);
    auto zMoveUL2 = new TGeoCombiTrans(0, 0, constants::tile::width * 2., rot);
    auto zMoveULDB = new TGeoCombiTrans(0, 0, -databackbone::width / 2., rot);
    // Lets attach the tiles to the QS.
    mRSU->AddNode(mTile, nCopyRSU++, rot);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveUL1);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveUL2);
    mRSU->AddNode(dataBackboneVol, nCopyDB++, zMoveULDB);

    // Upper Right
    auto zMoveUR0 = new TGeoCombiTrans(0, 0, +width / 2., rot);
    auto zMoveUR1 = new TGeoCombiTrans(0, 0, constants::tile::width + width / 2., rot);
    auto zMoveUR2 = new TGeoCombiTrans(0, 0, constants::tile::width * 2. + width / 2., rot);
    auto zMoveURDB = new TGeoCombiTrans(0, 0, -databackbone::width / 2. + width / 2., rot);
    // Lets attach the tiles to the QS.
    mRSU->AddNode(mTile, nCopyRSU++, zMoveUR0);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveUR1);
    mRSU->AddNode(mTile, nCopyRSU++, zMoveUR2);
    mRSU->AddNode(dataBackboneVol, nCopyDB++, zMoveURDB);
  }

  void createSegment()
  {
    using namespace constants::segment;
    // A segment is 12 RSUs + left and right end cap. We place the first rsu
    // as z-coordinate center and attach to this. Hence, we will displace the
    // left end-cap to the left and the right to right.
    mSegment = new TGeoVolumeAssembly(Form("segment_%d", mNLayer));
    mSegment->VisibleDaughters();

    for (unsigned int i{0}; i < nRSUs; ++i) {
      auto zMove = new TGeoTranslation(0, 0, +i * constants::rsu::width + constants::rsu::databackbone::width);
      mSegment->AddNode(mRSU, i, zMove);
    }

    // LEC
    double lecPhi1 = 0;
    double lecPhi2 = lec::length / mRmin * constants::Rad2Deg;
    auto zMoveLEC = new TGeoTranslation(0, 0, -lec::width / 2.);
    auto lec =
      new TGeoTubeSeg(mRmin, mRmax, lec::width / 2., lecPhi1, lecPhi2);
    auto lecVol = new TGeoVolume(
      Form("lec_%d", mNLayer), lec,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    lecVol->SetLineColor(lec::color);
    lecVol->RegisterYourself();
    if (mVerbose) {
      std::cout << "LEC:" << std::endl;
      lecVol->InspectShape();
      lecVol->InspectMaterial();
    }
    mSegment->AddNode(lecVol, 0, zMoveLEC);

    // REC; reuses lecPhi1,2
    auto zMoveREC = new TGeoTranslation(0, 0, nRSUs * constants::rsu::width + rec::width / 2.);
    auto rec =
      new TGeoTubeSeg(mRmin, mRmax, rec::width / 2., lecPhi1, lecPhi2);
    auto recVol = new TGeoVolume(
      Form("rec_%d", mNLayer), rec,
      gGeoManager->GetMedium(material::MaterialNames[material::DeadZone]));
    recVol->SetLineColor(rec::color);
    recVol->RegisterYourself();
    if (mVerbose) {
      std::cout << "REC:" << std::endl;
      recVol->InspectShape();
      recVol->InspectMaterial();
    }
    mSegment->AddNode(recVol, 0, zMoveREC);
  }

  void createChip()
  {

    // A HalfLayer is composed out of multiple segment stitched together along
    // rphi.
    mChip = new TGeoVolumeAssembly(Form("chip_%d", mNLayer));
    mChip->VisibleDaughters();

    for (unsigned int i{0}; i < constants::nSegments[mNLayer]; ++i) {
      double phiOffset = constants::segment::length / mRmin * constants::Rad2Deg;
      auto rot = new TGeoRotation("", 0, 0, phiOffset * i);
      mChip->AddNode(mSegment, i, rot);
    }
  }

  void createCarbonForm()
  {
    mCarbonForm = new TGeoVolumeAssembly(Form("carbonform_%d", mNLayer));
    mCarbonForm->SetVisibility(true);
    mCarbonForm->VisibleDaughters();
    using namespace constants::carbonfoam;
    double dRadius = -1;
    if (mNLayer < 2)
      dRadius = constants::radii[mNLayer + 1] - constants::radii[mNLayer] - constants::thickness;
    else
      dRadius = 0.7; // TODO: lack of carbon foam radius for layer 2, use 0.7mm as a temporary value

    double phiSta = edgeBetwChipAndFoam / (0.5 * constants::radii[mNLayer + 1] + constants::radii[mNLayer]) * constants::Rad2Deg;
    double phiEnd = (constants::nSegments[mNLayer] * constants::segment::length) / constants::radii[mNLayer] * constants::Rad2Deg - phiSta;
    double phiLongeronsCover = longeronsHeight / (0.5 * constants::radii[mNLayer + 1] + constants::radii[mNLayer]) * constants::Rad2Deg;

    // H-rings foam
    auto HringC = new TGeoTubeSeg(Form("HringC_%d", mNLayer), mRmax, mRmax + dRadius, HringWidth / 2., phiSta, phiEnd);
    auto HringA = new TGeoTubeSeg(Form("HringA_%d", mNLayer), mRmax, mRmax + dRadius, HringWidth / 2., phiSta, phiEnd);
    auto HringCWithHoles = getHringShape(HringC, nHoles[mNLayer], radiusHoles[mNLayer]);
    auto HringAWithHoles = getHringShape(HringA, nHoles[mNLayer], radiusHoles[mNLayer]);
    auto HringCVol = new TGeoVolume(Form("hringC_%d", mNLayer), HringCWithHoles, gGeoManager->GetMedium(material::MaterialNames[material::Silicon]));
    HringCVol->SetLineColor(color);
    auto HringAVol = new TGeoVolume(Form("hringA_%d", mNLayer), HringAWithHoles, gGeoManager->GetMedium(material::MaterialNames[material::Silicon]));
    HringAVol->SetLineColor(color);
    auto zMoveHringC = new TGeoTranslation(0, 0, -constants::segment::lec::width + HringWidth / 2.);
    auto zMoveHringA = new TGeoTranslation(0, 0, -constants::segment::lec::width + HringWidth / 2. + constants::segment::width - HringWidth);

    // Longerons are made by same material
    auto longeronR = new TGeoTubeSeg(Form("longeronR_%d", mNLayer), mRmax, mRmax + dRadius, longeronsWidth / 2, phiSta, phiSta + phiLongeronsCover);
    auto longeronL = new TGeoTubeSeg(Form("longeronL_%d", mNLayer), mRmax, mRmax + dRadius, longeronsWidth / 2, phiEnd - phiLongeronsCover, phiEnd);
    TString nameLongerons = Form("longeronR_%d + longeronL_%d", mNLayer, mNLayer);
    auto longerons = new TGeoCompositeShape(nameLongerons);
    auto longeronsVol = new TGeoVolume(Form("longerons_%d", mNLayer), longerons, gGeoManager->GetMedium(material::MaterialNames[material::Silicon]));
    longeronsVol->SetLineColor(color);
    auto zMoveLongerons = new TGeoTranslation(0, 0, -constants::segment::lec::width + constants::segment::width / 2.);

    mCarbonForm->AddNode(mChip, 0);
    mCarbonForm->AddNode(HringCVol, 0, zMoveHringC);
    mCarbonForm->AddNode(HringAVol, 0, zMoveHringA);
    mCarbonForm->AddNode(longeronsVol, 0, zMoveLongerons);
  }

  TGeoCompositeShape* getHringShape(TGeoTubeSeg* Hring, int nHoles, double radiusHoles)
  {
    //////////////////////////////////////
    // Function to dig holes in H-rings //
    //////////////////////////////////////
    double stepPhiHoles = (Hring->GetPhi2() - Hring->GetPhi1()) / (nHoles);
    double phiHolesSta = Hring->GetPhi1() + stepPhiHoles / 2.;
    double radiusHring = 0.5 * (Hring->GetRmin() + Hring->GetRmax());
    TGeoCompositeShape* HringWithHoles = nullptr;
    TString nameAllHoles = "";
    for (int iHoles = 0; iHoles < nHoles; iHoles++) {
      double phiHole = phiHolesSta + stepPhiHoles * iHoles;
      TString nameHole = Form("hole_%d_%d", iHoles, mNLayer);
      auto hole = new TGeoTube(nameHole, 0, radiusHoles, 3 * Hring->GetDz());
      // move hole to the hring radius
      auto zMoveHole = new TGeoTranslation(Form("zMoveHole_%d_%d", iHoles, mNLayer), radiusHring * cos(phiHole * constants::Deg2Rad), radiusHring * sin(phiHole * constants::Deg2Rad), 0);
      zMoveHole->RegisterYourself();
      nameAllHoles += Form("hole_%d_%d:zMoveHole_%d_%d + ", iHoles, mNLayer, iHoles, mNLayer);
    }
    nameAllHoles.Remove(nameAllHoles.Length() - 3, 3);
    TString nameHringWithHoles = Form("%s - (%s)", Hring->GetName(), nameAllHoles.Data());
    HringWithHoles = new TGeoCompositeShape(nameHringWithHoles);
    return HringWithHoles;
  }
  void createLayerImpl()
  {
    // At long last a single layer... A layer is two HalfLayers (duuhhh) but
    // we have to take care of the equatorial gap. So both half layers will be
    // offset slightly by rotating in phi the upper HalfLayer and negative phi
    // the other one.
    mLayer = new TGeoVolumeAssembly(Form("layer_%d", mNLayer));
    mLayer->VisibleDaughters();

    // The offset is the right angle triangle of the middle radius with the
    // transverse axis.
    double phiOffset = std::asin(constants::equatorialGap / mR) * constants::Rad2Deg;
    // double phiOffset = constants::equatorialGap / mRmin / 2.;
    auto rotTop = new TGeoRotation("", 0, 0, +phiOffset);
    auto rotBot = new TGeoRotation("", 0, 0, phiOffset + 180);

    mLayer->AddNode(mCarbonForm, 0, rotTop);
    mLayer->AddNode(mCarbonForm, 1, rotBot);
  }

 private:
  bool mVerbose{false};   // Verbose debug output
  bool mSubstrate{false}; // create substrate layer
  uint8_t mNLayer{0};     // Layer number
  double mRmin{0};        // Minimum Radius
  double mR{0};           // Middle Radius
  double mRmax{0};        // Maximum Radius

  // Individual Pieces
  TGeoVolume* mPixelArray{nullptr};
  TGeoVolumeAssembly* mTile{nullptr};
  TGeoVolumeAssembly* mRSU{nullptr};
  TGeoVolumeAssembly* mSegment{nullptr};
  TGeoVolumeAssembly* mChip{nullptr};
  TGeoVolumeAssembly* mCarbonForm{nullptr};
  TGeoVolumeAssembly* mLayer{nullptr};

  ClassDefNV(ITS3Layer, 0);
};

template <typename value_t = double>
class SegmentationSuperAlpide
{
  // This class defines the segmenation of the pixelArray in the tile. We define
  // two coordinate systems, one width x,z detector local coordianates (cm) and
  // the more natural row,col layout: Also all the transformation between these
  // two. The class provides the transformation from the tile to TGeo
  // coordinates.

  // row,col=0
  // |
  // v
  // x----------------------x
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |                        ^ x
  // |           |          |                        |
  // |           |          |                        |
  // |           |          |                        |
  // |-----------X----------|  X marks (x,z)=(0,0)   X----> z
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // x----------------------x

 public:
  SegmentationSuperAlpide(int layer = 0) : mLayer{layer} {}

  /// Transformation from the curved surface to a flat surface
  /// It works only if the detector is not rototraslated
  /// \param xCurved Detector local curved coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yCurved Detector local curved coordinate y in cm with respect to
  /// the center of the sensitive volume.
  /// \param xFlat Detector local flat coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yFlat Detector local flat coordinate y in cm with respect to
  /// the center of the sensitive volume.
  void curvedToFlat(value_t xCurved, value_t yCurved, value_t& xFlat, value_t& yFlat)
  {
    // FIXME: Chunzheng: I think we should use the phi instead of complementary angles of phi
    // value_t dist = std::hypot(xCurved, yCurved);
    // yFlat = dist - mEffRadius;
    // value_t phi = (value_t)constants::PI / 2 - std::atan2((double)yCurved, (double)xCurved);
    // xFlat = mEffRadius * phi;

    value_t dist = std::hypot(xCurved, yCurved);
    yFlat = dist - mEffRadius;
    value_t phi = std::acos(xCurved / dist);
    xFlat = dist * phi;
  }

  /// Transformation from the flat surface to a curved surface
  /// It works only if the detector is not rototraslated
  /// \param xFlat Detector local flat coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yFlat Detector local flat coordinate y in cm with respect to
  /// the center of the sensitive volume.
  /// \param xCurved Detector local curved coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yCurved Detector local curved coordinate y in cm with respect to
  /// the center of the sensitive volume.
  void flatToCurved(value_t xFlat, value_t yFlat, value_t& xCurved, value_t& yCurved)
  {
    // FIXME: Chunzheng: I think maybe the xCurved and yCurved should be swapped in the old code, need to check
    // value_t dist = yFlat + mEffRadius;
    // value_t phi = xFlat / dist;
    // value_t tang = std::tan((value_t)constants::PI / 2 - (value_t)phi);
    // xCurved = (xFlat > 0 ? 1.f : -1.f) * dist / std::sqrt(1 + tang * tang);
    // yCurved = xCurved * tang;

    // we don't need to calculate the tangent
    // when flat to curved, xflat is the arc length, so we can just get the Central angle easily
    value_t dist = yFlat + mEffRadius;
    value_t phi = xFlat / dist;
    xCurved = dist * std::cos(phi);
    yCurved = dist * std::sin(phi);
  }

  /// Transformation from Geant detector centered local coordinates (cm) to
  /// Pixel cell numbers iRow and iCol.
  /// Returns true if point x,z is inside sensitive volume, false otherwise.
  /// A value of -1 for iRow or iCol indicates that this point is outside of the
  /// detector segmentation as defined.
  /// \param float x Detector local coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param float z Detector local coordinate z in cm with respect to
  /// the center of the sensitive volume.
  /// \param int iRow Detector x cell coordinate.
  /// \param int iCol Detector z cell coordinate.
  bool localToDetector(value_t const xRow, value_t const zCol, int& iRow, int& iCol) const noexcept
  {
    localToDetectorUnchecked(xRow, zCol, iRow, iCol);
    if (!isValid(iRow, iCol)) {
      iRow = iCol = -1;
      return false;
    }
    return true;
  }

  // Same as localToDetector w.o. checks.
  void localToDetectorUnchecked(value_t const xRow, value_t const zCol, int& iRow, int& iCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    // value_t x = cp::length / 2. - xRow; // transformation to upper edge of pixelarray
    // value_t z = zCol + cp::width / 2.;  // transformation to left edge of pixelarray
    // iRow = std::floor(x / cp::pixel::pitchRow);
    // iCol = std::floor(z / cp::pixel::pitchCol);

    // @Chunzheng
    iRow = std::floor((cp::length / 2. - xRow) / cp::pixel::pitchRow);
    iCol = std::floor((zCol + cp::width / 2.) / cp::pixel::pitchCol);
  }

  /// Transformation from Detector cell coordinates to Geant detector centered
  /// local coordinates (cm)
  /// \param int iRow Detector x cell coordinate.
  /// \param int iCol Detector z cell coordinate.
  /// \param float x Detector local coordinate x in cm with respect to the
  /// center of the sensitive volume.
  /// \param float z Detector local coordinate z in cm with respect to the
  /// center of the sensitive volume.
  /// If iRow and or iCol is outside of the segmentation range a value of -0.5*Dx()
  /// or -0.5*Dz() is returned.
  bool detectorToLocal(int const iRow, int const iCol, value_t& xRow, value_t& zCol) const noexcept
  {
    detectorToLocalUnchecked(iRow, iCol, xRow, zCol);
    if (!isValid(xRow, zCol)) {
      return false;
    }
    return true;
  }

  // Same as detectorToLocal w.o. checks.
  // We position ourself in the middle of the pixel.
  void detectorToLocalUnchecked(int iRow, int iCol, value_t& xRow, value_t& zCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    // xRow = -(iRow - 0.5) * cp::pixel::pitchRow + cp::length / 2.;
    // zCol = -(iCol + 0.5) * cp::pixel::pitchCol - cp::width / 2.;

    // @Chunzheng
    xRow = -(iRow + 0.5) * cp::pixel::pitchRow + cp::length / 2.;
    zCol = (iCol + 0.5) * cp::pixel::pitchCol - cp::width / 2.;
  }

  bool detectorToLocal(float row, float col, float& xRow, float& zCol)
  {
    return detectorToLocal(static_cast<int>(row), static_cast<int>(col), xRow, zCol);
  }

  void detectorToLocalUnchecked(float row, float col, float& xRow, float& zCol)
  {
    detectorToLocalUnchecked(static_cast<int>(row), static_cast<int>(col), xRow, zCol);
  }

 private:
  bool isValid(value_t const xRow, value_t const zCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    if (xRow < 0. || xRow >= cp::length || zCol < 0. || zCol >= cp::width) {
      return false;
    }
    return true;
  }

  bool isValid(int const iRow, int const iCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    if (iRow < 0 || iRow >= cp::nRows || iCol < 0 || iCol >= cp::nCols) {
      return false;
    }
    return true;
  }

  const int mLayer; ///< chip layer
  const value_t mEffRadius{constants::radii[mLayer] + constants::thickness / 2.};
}; // namespace its3

} // namespace its3
} // namespace my

void tileSegTest()
{
  using namespace my::its3;
  gSystem->Load("libGeom");
  new TGeoManager("geo", "poza7");
  material::Material::Instance().print();
  TGeoVolume* top = gGeoManager->MakeBox(
    "TOP", gGeoManager->GetMedium(material::MaterialNames[material::Vacuum]),
    100., 100., 100.);
  gGeoManager->SetTopVolume(top);

  ITS3Layer layer;
  layer.createLayer(top, 0, true);
  layer.createLayer(top, 1, true);
  layer.createLayer(top, 2, true);
  gGeoManager->CloseGeometry();
  gGeoManager->SetVisLevel(10);
  /* gGeoManager->SetMaxVisNodes(20000); */
  /* gGeoManager->SetVisOption(1); */
  gGeoManager->OptimizeVoxels();
  /* gGeoManager->CheckGeometryFull(); */
  /* gGeoManager->RandomRays(); */

  TFile* file = TFile::Open("tile_seg.root", "RECREATE");
  gGeoManager->Write();
  file->Close();

  // Event Display
  if (false) {
    TEveManager::Create();
    TGeoNode* tnode = gGeoManager->GetTopNode();
    tnode->SetVisibility(kFALSE);
    TEveGeoTopNode* eve_tnode = new TEveGeoTopNode(gGeoManager, tnode);
    eve_tnode->SetVisLevel(0);
    gEve->AddGlobalElement(eve_tnode);

    auto list = new TEveTrackList();
    list->SetName("Heix Propagator");
    auto prop = list->GetPropagator();
    prop->SetFitDaughters(kFALSE);
    prop->SetMaxZ(20);
    prop->SetMaxR(5);
    prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., -0.5));
    list->SetElementName(Form("%s, constB", list->GetElementName()));
    auto rc = new TEveRecTrackD();
    rc->fV.Set(0.028558, -0.000918, 3.691919);
    rc->fP.Set(0.767095, -2.400006, -0.313103);
    rc->fSign = -1;
    auto track = new TEveTrack(rc, prop);
    track->SetName(Form("Charge %d", -1));
    list->SetLineColor(kCyan);
    track->SetLineColor(list->GetLineColor());

    gEve->AddElement(list);
    list->AddElement(track);

    track->MakeTrack();

    TEveViewer* ev = gEve->GetDefaultViewer();
    TGLViewer* gv = ev->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, true, true, nullptr);
    gEve->FullRedraw3D(kTRUE);
    gSystem->ProcessEvents();

    gv->CurrentCamera().RotateRad(-0.5, 1.4);
    gv->RequestDraw();
  } else {
    gGeoManager->GetTopVolume()->Draw("ogl");
  }

  if (true) {
    namespace cp = constants::pixelarray;
    TH2I* h_raw_col = new TH2I("h_raw_col", "h_raw_col;raw;col", cp::nRows, 0, cp::nRows, cp::nCols, 0, cp::nCols);
    TH2D* h_xLocal_zLocal = new TH2D("h_xLocal_zLocal", "h_xLocal_zLocal;xLocal;yLocal", 1000, -cp::length / 2, cp::length / 2, 1000, -cp::width / 2, cp::width / 2);
    TH2D* h_xCurved_yCurved = new TH2D("h_xCurved_yCurved", "h_xCurved_yCurved", 100, -cp::length / 2, cp::length / 2, 100, -cp::width / 2, cp::width / 2);
    TGraph* g_raw_xLocal = new TGraph();
    TGraph* g_col_zLocal = new TGraph();

    TH2D* h_xGlobal_yGlobal = new TH2D("h_xGlobal_yGlobal", "h_xGlobal_yGlobal;xGlobel;yGlobal", 1000, -constants::radii[2] * 1.5, +constants::radii[2] * 1.5, 1000, -constants::radii[2] * 1.5, +constants::radii[2] * 1.5);
    TH2D* h_zGlobal_xGlobal = new TH2D("h_zGlobal_xGlobal", "h_zGlobal_xGlobal;xGlobel;yGlobal", 1000, -constants::segment::width, +constants::segment::width, 100, -constants::radii[2] * 1.5, +constants::radii[2] * 1.5);

    for (unsigned int iLayer{}; iLayer < constants::nLayers; ++iLayer) {
      for (unsigned int iCarbonForm{0}; iCarbonForm < 2; ++iCarbonForm) {
        // No Loop for chip = carbonform id
        for (unsigned int iSegment{0}; iSegment < constants::nSegments[iLayer]; ++iSegment) {
          for (unsigned int iRSU{0}; iRSU < 12; ++iRSU) {
            for (unsigned int iTile{0}; iTile < 12; ++iTile) {
              // auto path = fmt::format("/TOP_1/layer_{}_0/carbonform_{}_{}/chip_{}_0/segment_{}_{}/rsu_{}_{}/tile_{}_{}/pixelarray_{}_0",
              //             iLayer, iLayer, iCarbonForm, iLayer, iLayer, iSegment, iLayer, iRSU, iLayer, iTile, iLayer);
              TString path = Form("/TOP_1/layer_%d_0/carbonform_%d_%d/chip_%d_0/segment_%d_%d/rsu_%d_%d/tile_%d_%d/pixelarray_%d_0",
                                  iLayer, iLayer, iCarbonForm, iLayer, iLayer, iSegment, iLayer, iRSU, iLayer, iTile, iLayer);
              if (!gGeoManager->CheckPath(path.Data())) {
                std::cerr << path << std::endl;
              }
              gGeoManager->cd(path.Data());
              auto pixelArray = gGeoManager->GetCurrentVolume();
              // Get the current matrix
              TGeoHMatrix* matrix = gGeoManager->GetCurrentMatrix();
              int iRow = 0;
              int iCol = 0;
              double xLocal = 0;
              double zLocal = 0;
              double xCurved = 0;
              double yCurved = 0;

              // Test the coordinate from dectector(row,col) to local(x',z') to curved(x'',y'') to global(x,y,z)
              for (size_t i = 0; i < 1000; i++) {
                //randomly sow the points in the pixel array
                int row = gRandom->Uniform(0, cp::nRows);
                int col = gRandom->Uniform(0, cp::nCols);
                double xLocal = 0;
                double zLocal = 0;
                double xCurved = 0;
                double yCurved = 0;

                SegmentationSuperAlpide<double> seg(iLayer);
                seg.detectorToLocal(row, col, xLocal, zLocal);
                seg.flatToCurved(xLocal, 0, xCurved, yCurved);
                double posLocal[3] = {xCurved, yCurved, zLocal};
                double posGlobal[3] = {0, 0, 0};
                matrix->LocalToMaster(posLocal, posGlobal);

                h_raw_col->Fill(row, col);
                h_xLocal_zLocal->Fill(xLocal, zLocal);
                g_raw_xLocal->SetPoint(i, row, xLocal);
                g_col_zLocal->SetPoint(i, col, zLocal);
                h_xCurved_yCurved->Fill(xCurved, yCurved);
                h_xGlobal_yGlobal->Fill(posGlobal[0], posGlobal[1]);
                h_zGlobal_xGlobal->Fill(posGlobal[2], posGlobal[0]);
              }
            }
          }
        }
      }
    }
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    gStyle->SetPadLeftMargin(0.15);
    c1->Divide(2, 2);
    c1->cd(1);
    h_raw_col->Draw("colz");
    c1->cd(2);
    h_xLocal_zLocal->Draw("colz");
    c1->cd(3);
    g_raw_xLocal->SetTitle(";raw;xLocal");
    g_raw_xLocal->SetMarkerStyle(20);
    g_raw_xLocal->SetMarkerSize(0.2);
    g_raw_xLocal->Draw("ap");
    c1->cd(4);
    g_col_zLocal->SetTitle(";col;zLocal");
    g_col_zLocal->Draw("ap");

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
    h_xCurved_yCurved->Draw("colz");

    TCanvas* c3 = new TCanvas("c3", "c3", 800, 400);
    c3->Divide(2, 1);
    c3->cd(1);
    h_xGlobal_yGlobal->Draw("colz");
    c3->cd(2);
    h_zGlobal_xGlobal->Draw("colz");
  }
}
