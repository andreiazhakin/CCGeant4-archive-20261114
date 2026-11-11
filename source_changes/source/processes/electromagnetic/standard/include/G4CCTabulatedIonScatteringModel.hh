#ifndef G4CC_TABULATED_ION_SCATTERING_MODEL_HH
#define G4CC_TABULATED_ION_SCATTERING_MODEL_HH

#include "G4VEmModel.hh"
#include "G4CCElasticTableData.hh"
#include "G4ThreeVector.hh"
#include <cstdint>
#include <memory>

class G4Track;

namespace tabulated {


struct ElasticKinematics {
  double finalKineticEnergy = 0;
  double recoilKineticEnergy = 0;
  double momentumCM = 0;
  double q2 = 0;
  G4ThreeVector primaryMomentum;
  G4ThreeVector recoilMomentum;
};

ElasticKinematics ComputeElasticKinematics(double projectileMass,
                                           double targetMass,
                                           double kineticEnergy,
                                           double x,
                                           double phi);

struct CollisionRecord {
  std::uint64_t sequence = 0;
  int trackId = -1;
  int projectileZ = 0;
  int projectileA = 0;
  int targetZ = 0;
  int targetA = 0;
  int targetIsomer = 0;
  int nodeIndex = -1;
  bool tabulated = false;
  bool nativeDelegated = false;
  bool recoilCreated = false;
  double incidentT_MeV = 0;
  double incidentCharge_e = 0;
  double x = 0;
  double q2_MeV2 = 0;
  double thetaCM_rad = 0;
  double recoilT_MeV = 0;
  double finalT_MeV = 0;
  double energyResidual_MeV = 0;
  double momentumResidual_MeV = 0;
  double targetMass_MeV = 0;
};

class G4CCTabulatedIonScatteringModel final : public G4VEmModel {
 public:
  explicit G4CCTabulatedIonScatteringModel(
      std::shared_ptr<const ElasticTableData> data, bool useTables = true);
  ~G4CCTabulatedIonScatteringModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  void InitialiseLocal(const G4ParticleDefinition*, G4VEmModel*) override;
  void StartTracking(G4Track*) override;

  G4double CrossSectionPerVolume(const G4Material*, const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy = 0,
                                G4double maxEnergy = DBL_MAX) override;
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                     G4double kineticEnergy, G4double Z,
                                     G4double A = 0,
                                     G4double cutEnergy = 0,
                                     G4double maxEnergy = DBL_MAX) override;
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
                         const G4DynamicParticle*, G4double = 0,
                         G4double = DBL_MAX) override;
  void ModelDescription(std::ostream&) const override;

  void SetRecoilThreshold(G4double);
  G4double CrossSectionForIsotope(const G4ParticleDefinition*, G4double energy,
                                 int Z, int A, int isomerLevel = 0);
  static const CollisionRecord& LastCollision();
  static void ResetThreadDiagnostics();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

using TabulatedIonElasticModel = G4CCTabulatedIonScatteringModel;
}

using G4CCTabulatedIonScatteringModel =
    tabulated::G4CCTabulatedIonScatteringModel;
#endif
