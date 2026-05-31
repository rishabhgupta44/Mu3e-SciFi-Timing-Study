// Copyright (c) Rishabh Gupta 2026
// This is a test code used for learning.

#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandPoisson.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>

//this config is based on arXiv:2501.14692 
namespace Config {
constexpr G4double fiberDiameter = 250.0 * micrometer;
constexpr G4double fiberRadius = 0.5 * fiberDiameter;
constexpr G4double fiberLength = 100.0 * mm;

constexpr G4int nLayers = 3;
constexpr G4int nFibersPerLayer = 31;

constexpr G4double pitch = fiberDiameter;
constexpr G4double layerSpacing = 0.5 * std::sqrt(3.0) * fiberDiameter;

constexpr G4double worldX = 12.0 * mm;
constexpr G4double worldY = 12.0 * mm;
constexpr G4double worldZ = 70.0 * mm;

constexpr G4double lightYield = 8000.0;
constexpr G4double trappingEfficiency = 0.055;
constexpr G4double attenuationLength = 1500.0 * mm;
constexpr G4double photonDetectionEfficiency = 0.25;
constexpr G4double refractiveIndexEff = 1.59;
constexpr G4double sigmaPhotonTerm = 0.80 * ns;
constexpr G4double sigmaElectronics = 0.07 * ns;
constexpr G4double timeWalkCoeff = 0.20 * ns;
constexpr G4double thresholdNpe = 3.0;
constexpr G4double birksK = 0.126 * mm / MeV;
}

struct HitAccum {
  G4double evis = 0.0;
  G4double edep = 0.0;
  G4ThreeVector weightedPos = G4ThreeVector();
  G4double weightedTime = 0.0;
};

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  G4LogicalVolume* fiberLogic = nullptr;

  G4VPhysicalVolume* Construct() override {
    auto nist = G4NistManager::Instance();
    auto air = nist->FindOrBuildMaterial("G4_AIR");
    auto polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    auto solidWorld = new G4Box("World", Config::worldX, Config::worldY, Config::worldZ);
    auto logicWorld = new G4LogicalVolume(solidWorld, air, "WorldLV");

    auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "WorldPV", nullptr, false, 0, true);

    auto solidFiber = new G4Tubs(
        "SciFiFiberSolid", 0.0, Config::fiberRadius, 0.5 * Config::fiberLength, 0.0, 360.0 * deg);

    const G4Colour layerColours[] = {
        G4Colour(0.85, 0.20, 0.20, 0.75),
        G4Colour(0.20, 0.70, 0.25, 0.75),
        G4Colour(0.20, 0.35, 0.90, 0.75)};

    G4LogicalVolume* layerLogic[Config::nLayers] = {};
    for (G4int layer = 0; layer < Config::nLayers; ++layer) {
      layerLogic[layer] = new G4LogicalVolume(solidFiber, polystyrene, "SciFiFiberLV");
      auto* vis = new G4VisAttributes(layerColours[layer]);
      vis->SetForceSolid(true);
      layerLogic[layer]->SetVisAttributes(vis);
    }

    fiberLogic = layerLogic[0];

    const G4double firstX = -0.5 * (Config::nFibersPerLayer - 1) * Config::pitch;

    for (G4int layer = 0; layer < Config::nLayers; ++layer) {
      const G4double y = (layer - 1) * Config::layerSpacing;
      const G4double stagger = (layer == 1) ? 0.5 * Config::pitch : 0.0;

      for (G4int i = 0; i < Config::nFibersPerLayer; ++i) {
        const G4double x = firstX + i * Config::pitch + stagger;
        const G4int copyNo = layer * 1000 + i;

        new G4PVPlacement(
            nullptr,
            G4ThreeVector(x, y, 0.0),
          layerLogic[layer],
            "SciFiFiberPV",
            logicWorld,
            false,
            copyNo,
            true);
      }
    }

    return physWorld;
  }
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 private:
  G4ParticleGun* gun = nullptr;

 public:
  PrimaryGeneratorAction() {
    gun = new G4ParticleGun(1);
    auto particle = G4ParticleTable::GetParticleTable()->FindParticle("e+");
    gun->SetParticleDefinition(particle);
    gun->SetParticleMomentumDirection(G4ThreeVector(0.0, 1.0, 0.0));
  }

  ~PrimaryGeneratorAction() override { delete gun; }

  void GeneratePrimaries(G4Event* event) override {
    const G4double energy = (10.0 + 43.0 * G4UniformRand()) * MeV;

    const G4double ribbonHalfWidth = 0.5 * Config::nFibersPerLayer * Config::pitch;
    const G4double x = (2.0 * G4UniformRand() - 1.0) * ribbonHalfWidth;
    const G4double z = (2.0 * G4UniformRand() - 1.0) * 0.45 * Config::fiberLength;
    const G4double y = -3.0 * mm;

    const G4double thetaX = (2.0 * G4UniformRand() - 1.0) * 20.0 * mrad;
    const G4double thetaZ = (2.0 * G4UniformRand() - 1.0) * 20.0 * mrad;

    G4ThreeVector dir(thetaX, 1.0, thetaZ);
    dir = dir.unit();

    gun->SetParticleEnergy(energy);
    gun->SetParticlePosition(G4ThreeVector(x, y, z));
    gun->SetParticleMomentumDirection(dir);
    gun->GeneratePrimaryVertex(event);
  }
};

class RunAction : public G4UserRunAction {
 public:
  static std::ofstream out;

  void BeginOfRunAction(const G4Run*) override {
    std::ifstream existing("scifi_hits.csv", std::ios::binary | std::ios::ate);
    const bool needsHeader = !existing.good() || existing.tellg() == 0;

    out.open("scifi_hits.csv", std::ios::app);
    if (needsHeader) {
      out << "eventID,layerID,fiberID,copyNo,"
          << "edep_keV,evis_keV,"
          << "x_um,y_um,z_mm,t_true_ns,"
          << "npe_left,npe_right,npe_total,"
          << "t_left_ns,t_right_ns,t_reco_ns,z_reco_mm,"
          << "dt_ps,dz_um,valid\n";
    }
  }

  void EndOfRunAction(const G4Run*) override { out.close(); }
};

std::ofstream RunAction::out;

class EventAction : public G4UserEventAction {
 private:
  std::map<G4int, HitAccum> hits;

 public:
  void BeginOfEventAction(const G4Event*) override { hits.clear(); }

  void AddHit(G4int copyNo, G4double edep, G4double evis, const G4ThreeVector& pos, G4double time) {
    auto& hit = hits[copyNo];
    hit.edep += edep;
    hit.evis += evis;
    hit.weightedPos += evis * pos;
    hit.weightedTime += evis * time;
  }

  void EndOfEventAction(const G4Event* event) override {
    const G4int eventID = event->GetEventID();
    const G4double vEff = c_light / Config::refractiveIndexEff;

    for (const auto& item : hits) {
      const G4int copyNo = item.first;
      const HitAccum& hit = item.second;

      if (hit.evis <= 0.0) continue;

      const G4int layerID = copyNo / 1000;
      const G4int fiberID = copyNo % 1000;

      const G4ThreeVector pos = hit.weightedPos / hit.evis;
      const G4double tTrue = hit.weightedTime / hit.evis;

      const G4double zFromLeft = pos.z() + 0.5 * Config::fiberLength;
      const G4double zToRight = 0.5 * Config::fiberLength - pos.z();

      const G4double nGammaMean = Config::lightYield * (hit.evis / MeV);
      const G4double nPhotLeftMean =
          nGammaMean * Config::trappingEfficiency * std::exp(-zFromLeft / Config::attenuationLength);
      const G4double nPhotRightMean =
          nGammaMean * Config::trappingEfficiency * std::exp(-zToRight / Config::attenuationLength);

      const G4double npeLeftMean = nPhotLeftMean * Config::photonDetectionEfficiency;
      const G4double npeRightMean = nPhotRightMean * Config::photonDetectionEfficiency;

      const G4double npeLeft = CLHEP::RandPoisson::shoot(std::max(0.0, npeLeftMean));
      const G4double npeRight = CLHEP::RandPoisson::shoot(std::max(0.0, npeRightMean));
      const G4double npeTotal = npeLeft + npeRight;

      const G4bool valid = (npeLeft >= Config::thresholdNpe) && (npeRight >= Config::thresholdNpe);

      G4double tLeft = -999.0 * ns;
      G4double tRight = -999.0 * ns;
      G4double tReco = -999.0 * ns;
      G4double zReco = -999.0 * mm;
      G4double dt = -999.0 * ns;
      G4double dz = -999.0 * mm;

      if (valid) {
        const G4double sigmaLeft = std::sqrt(
            std::pow(Config::sigmaPhotonTerm / std::sqrt(npeLeft), 2) + std::pow(Config::sigmaElectronics, 2));
        const G4double sigmaRight = std::sqrt(
            std::pow(Config::sigmaPhotonTerm / std::sqrt(npeRight), 2) + std::pow(Config::sigmaElectronics, 2));

        const G4double timeWalkLeft = Config::timeWalkCoeff / std::sqrt(npeLeft);
        const G4double timeWalkRight = Config::timeWalkCoeff / std::sqrt(npeRight);

        const G4double tLeftRaw = tTrue + zFromLeft / vEff + timeWalkLeft + G4RandGauss::shoot(0.0, sigmaLeft);
        const G4double tRightRaw = tTrue + zToRight / vEff + timeWalkRight + G4RandGauss::shoot(0.0, sigmaRight);

        tLeft = tLeftRaw;
        tRight = tRightRaw;

        const G4double tLeftCorr = tLeftRaw - timeWalkLeft;
        const G4double tRightCorr = tRightRaw - timeWalkRight;

        tReco = 0.5 * (tLeftCorr + tRightCorr - Config::fiberLength / vEff);
        zReco = 0.5 * vEff * (tLeftCorr - tRightCorr);
        dt = tReco - tTrue;
        dz = zReco - pos.z();
      }

      if (RunAction::out.is_open()) {
        RunAction::out << eventID << "," << layerID << "," << fiberID << "," << copyNo << ","
                       << hit.edep / keV << "," << hit.evis / keV << ","
                       << pos.x() / micrometer << "," << pos.y() / micrometer << "," << pos.z() / mm << ","
                       << tTrue / ns << ","
                       << npeLeft << "," << npeRight << "," << npeTotal << ","
                       << tLeft / ns << "," << tRight / ns << "," << tReco / ns << "," << zReco / mm << ","
                       << dt / ps << "," << dz / micrometer << "," << (valid ? 1 : 0) << "\n";
      }
    }
  }
};

class SteppingAction : public G4UserSteppingAction {
 private:
  EventAction* eventAction = nullptr;

 public:
  explicit SteppingAction(EventAction* action) : eventAction(action) {}

  void UserSteppingAction(const G4Step* step) override {
    const auto prePoint = step->GetPreStepPoint();
    const auto volume = prePoint->GetTouchableHandle()->GetVolume();
    if (!volume) return;

    if (volume->GetName() != "SciFiFiberPV") return;

    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    const G4double stepLength = step->GetStepLength();
    G4double evis = edep;

    if (stepLength > 0.0) {
      const G4double dEdx = edep / stepLength;
      evis = edep / (1.0 + Config::birksK * dEdx);
    }

    const G4int copyNo = prePoint->GetTouchableHandle()->GetCopyNumber();
    const G4ThreeVector pos = 0.5 * (prePoint->GetPosition() + step->GetPostStepPoint()->GetPosition());
    const G4double time = prePoint->GetGlobalTime();

    eventAction->AddHit(copyNo, edep, evis, pos, time);
  }
};

class ActionInitialization : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new PrimaryGeneratorAction());

    auto* runAction = new RunAction();
    SetUserAction(runAction);

    auto* eventAction = new EventAction();
    SetUserAction(eventAction);

    SetUserAction(new SteppingAction(eventAction));
  }
};

int main(int argc, char** argv) {
  auto* runManager = new G4RunManager();

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new FTFP_BERT());
  runManager->SetUserInitialization(new ActionInitialization());

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  auto* visManager = new G4VisExecutive();
  visManager->Initialize();

  auto* uiManager = G4UImanager::GetUIpointer();

  if (ui) {
    uiManager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    uiManager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;
  return 0;
}
