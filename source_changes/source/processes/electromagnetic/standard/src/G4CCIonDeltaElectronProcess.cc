#include "G4CCIonDeltaElectronProcess.hh"

#include "G4CCIonChargeState.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCMottCorrection.hh"
#include "G4CCScreening.hh"
#include "G4CCDiracPhysics.hh"
#include "G4CCDiracTable.hh"
#include "G4Element.hh"
#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCAtimaEquilibriumDedxTable.hh"
#include "G4GenericIon.hh"
#include "G4Exception.hh"
#include "G4DeltaAngle.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4EmCorrections.hh"
#include "G4ForceCondition.hh"
#include "G4Material.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VEmAngularDistribution.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <string>


namespace
{
G4bool G4CCEnvFlagValue(const char* value, G4bool defaultValue)
{
  if(value == nullptr) { return defaultValue; }
  const std::string s(value);
  if(s == "0" || s == "false" || s == "FALSE" ||
     s == "off" || s == "OFF" || s == "no" || s == "NO") {
    return false;
  }
  if(s == "1" || s == "true" || s == "TRUE" ||
     s == "on" || s == "ON" || s == "yes" || s == "YES") {
    return true;
  }
  return defaultValue;
}

G4bool G4CCEnvFlag(const char* name, G4bool defaultValue)
{
  return G4CCEnvFlagValue(std::getenv(name), defaultValue);
}

G4bool G4CCUseMottDeltaElectronsRuntime()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const char* envValue = std::getenv("G4CC_USE_MOTT_DELTA");
  static const G4bool envIsSet = (envValue != nullptr);
  static const G4bool cachedEnv = G4CCEnvFlagValue(envValue, true);
  if(envIsSet) { return cachedEnv; }
  return G4CCCurrentChargeConfig::Instance().UseMottDeltaElectrons();
}

G4bool G4CCUseMottDeltaElectronsAtEnergy(G4double kineticEnergyMeVPerU)
{
  return G4CCUseMottDeltaElectronsRuntime()
    && G4CCCurrentChargeConfig::MottEnabledForEnergy(kineticEnergyMeVPerU);
}

G4bool G4CCMottDiagnosticsEnabled()
{
  static const G4bool enabled = G4CCEnvFlag("G4CC_MOTT_DIAGNOSTICS", false);
  return enabled;
}

G4bool G4CCFastMottIntegralEnabled()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) { return true; }
  static const G4bool enabled = G4CCEnvFlag("G4CC_FAST_MOTT_INTEGRAL", true);
  return enabled;
}

G4bool G4CCUseDeltaAngularGeneratorRuntime()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) { return true; }
  static const char* envValue = std::getenv("G4CC_USE_DELTA_ANGLE");
  static const char* envValueAlt = std::getenv("G4CC_USE_DELTA_ANGULAR_GENERATOR");
  static const G4bool envIsSet = (envValue != nullptr) || (envValueAlt != nullptr);
  static const G4bool enabled = G4CCEnvFlagValue((envValue != nullptr) ? envValue : envValueAlt, true);
  return envIsSet ? enabled : true;
}

G4bool G4CCUseDeltaFormFactorRuntime()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const G4bool enabled = G4CCEnvFlag("G4CC_USE_DELTA_FORM_FACTOR", true);
  return enabled;
}

G4double G4CCEnvDouble(const char* name, G4double defaultValue)
{
  const char* value = std::getenv(name);
  if(value == nullptr || value[0] == '\0') { return defaultValue; }
  try {
    return std::stod(value);
  } catch(...) {
    return defaultValue;
  }
}

G4double G4CCClamp01(G4double value)
{
  return std::max(0.0, std::min(1.0, value));
}

G4double G4CCDeltaFormFactorMinimumEnergyMeVPerU()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4CCCurrentChargeConfig::High0MeVPerU();
  static const G4double threshold =
    std::max(0.0, G4CCEnvDouble("G4CC_DELTA_FORM_FACTOR_MIN_ENERGY_MEV_U", 1500.0));
  return threshold;
}

G4double G4CCDeltaFormFactorWeight(G4double kineticEnergyMeVPerU)
{
  if(!G4CCUseDeltaFormFactorRuntime()) { return 0.0; }
  return (kineticEnergyMeVPerU >= G4CCDeltaFormFactorMinimumEnergyMeVPerU()) ? 1.0 : 0.0;
}

G4double G4CCProjectileFormFactorParameter(const G4ParticleDefinition* particle, G4int zProjectile)
{
  if(particle == nullptr) { return 0.0; }
  const G4double mass = particle->GetPDGMass();
  if(mass <= 0.0) { return 0.0; }

  const G4double spin = particle->GetPDGSpin();
  G4double x = 0.8426*CLHEP::GeV;
  if(spin == 0.0 && mass < CLHEP::GeV) {
    x = 0.736*CLHEP::GeV;
  } else if(zProjectile > 1) {
    x /= G4NistManager::Instance()->GetA27(zProjectile);
  }
  return (x > 0.0) ? 2.0*CLHEP::electron_mass_c2/(x*x) : 0.0;
}

G4double G4CCProjectileMagMoment2(const G4ParticleDefinition* particle, G4double mass)
{
  if(particle == nullptr || mass <= 0.0) { return 0.0; }
  const G4double aMag = 1.0/(0.5*CLHEP::eplus*CLHEP::hbar_Planck*CLHEP::c_squared);
  const G4double magmom = particle->GetPDGMagneticMoment()*mass*aMag;
  return magmom*magmom - 1.0;
}

G4double G4CCProjectileFormFactorAcceptance(G4double deltaEnergy,
                                            G4double beta2,
                                            G4double tMax,
                                            G4double totalEnergy2,
                                            G4double spin,
                                            G4double projectileMass,
                                            G4double projectileFormFactor,
                                            G4double projectileMagMoment2)
{
  const G4double x = projectileFormFactor*deltaEnergy;
  if(x <= 1.e-6 || deltaEnergy <= 0.0 || tMax <= 0.0) { return 1.0; }

  const G4double x1 = 1.0 + x;
  G4double grej = 1.0/(x1*x1);
  if(spin > 0.0 && projectileMass > 0.0 && totalEnergy2 > 0.0) {
    const G4double f1 = 0.5*deltaEnergy*deltaEnergy/totalEnergy2;
    const G4double f = G4CCMottCorrection::G4BornShape(beta2, deltaEnergy,
                                                        tMax, totalEnergy2, spin);
    if(f > 0.0) {
      const G4double x2 = 0.5*CLHEP::electron_mass_c2*deltaEnergy/(projectileMass*projectileMass);
      grej *= (1.0 + projectileMagMoment2*(x2 - f1/f)/(1.0 + x2));
    }
  }
  return G4CCClamp01(grej);
}

G4double G4CCFallbackFormFactorAcceptance(G4int Zp,
                                          G4int Qe,
                                          G4double deltaEnergy,
                                          G4double qEffective2,
                                          G4double globalWeight,
                                          G4double formFactorWeight,
                                          G4double formFactorAcceptance)
{
  const G4double ffw = G4CCClamp01(formFactorWeight);
  if(ffw <= 0.0) { return 1.0; }
  const G4double w = G4CCClamp01(globalWeight);
  if(w >= 1.0) { return 1.0; }

  const G4double z2Current = G4CCScreening::DeltaChargeSquare(Zp, Qe, deltaEnergy);
  const G4double z2Fallback = std::max(0.0, qEffective2);
  const G4double denom = w*z2Current + (1.0 - w)*z2Fallback;
  if(denom <= 0.0) { return 1.0; }

  const G4double ff = 1.0 + ffw*(G4CCClamp01(formFactorAcceptance) - 1.0);
  const G4double numer = w*z2Current + (1.0 - w)*z2Fallback*ff;
  return G4CCClamp01(numer/denom);
}

G4int G4CCSelectRandomTargetZ(const G4Material* material)
{
  if(material == nullptr) { return 0; }
  const G4int nElements = static_cast<G4int>(material->GetNumberOfElements());
  const auto* elements = material->GetElementVector();
  if(nElements <= 0 || elements == nullptr) { return 0; }
  if(nElements == 1) { return (*elements)[0]->GetZasInt(); }

  const G4double* atomDensities = material->GetVecNbOfAtomsPerVolume();
  if(atomDensities == nullptr) { return (*elements)[0]->GetZasInt(); }

  G4double totalElectronDensityWeight = 0.0;
  for(G4int i = 0; i < nElements; ++i) {
    totalElectronDensityWeight += atomDensities[i]*(*elements)[i]->GetZ();
  }
  if(totalElectronDensityWeight <= 0.0) { return (*elements)[0]->GetZasInt(); }

  const G4double sample = G4UniformRand()*totalElectronDensityWeight;
  G4double cumulative = 0.0;
  for(G4int i = 0; i < nElements; ++i) {
    cumulative += atomDensities[i]*(*elements)[i]->GetZ();
    if(sample <= cumulative) { return (*elements)[i]->GetZasInt(); }
  }
  return (*elements)[nElements - 1]->GetZasInt();
}

G4String G4CCDefaultCorrectionZMode()
{
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_CORRECTION_Z_MODE");
    if(value == nullptr) { value = std::getenv("G4CC_CORRECTION_CHARGE_MODE"); }
    if(value == nullptr) { return G4String("effective"); }
    return G4String(value);
  }();
  return mode;
}

G4String G4CCMottZMode()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4String("effective");
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_MOTT_Z_MODE");
    if(value == nullptr) { return G4CCDefaultCorrectionZMode(); }
    return G4String(value);
  }();
  return mode;
}

G4double G4CCEffectiveChargeInEplus(const G4ParticleDefinition* p,
                                    const G4Material* mat,
                                    G4double kineticEnergy)
{
  if(p == nullptr || mat == nullptr) { return 0.0; }
  auto* corrections = G4LossTableManager::Instance()->EmCorrections();
  if(corrections == nullptr) { return 0.0; }
  return std::abs(corrections->GetParticleCharge(p, mat, kineticEnergy)/CLHEP::eplus);
}

G4double G4CCEffectiveChargeSquareForDelta(const G4Track& track)
{
  G4CCEffectiveChargeInEplus(track.GetParticleDefinition(), track.GetMaterial(),
                             track.GetKineticEnergy());
  return G4CCIonChargeState::EffectiveChargeSquare(track);
}

G4int G4CCChooseMottZ(const G4Track& track, G4int zProjectile)
{
  if(zProjectile <= 0) { return 0; }

  const G4String mode = G4CCMottZMode();
  G4double z = static_cast<G4double>(zProjectile);

  if(mode == "effective" || mode == "eff" || mode == "qeff") {
    const G4double zeff = G4CCEffectiveChargeInEplus(track.GetParticleDefinition(),
                                                     track.GetMaterial(),
                                                     track.GetKineticEnergy());
    if(zeff > 0.0) { z = zeff; }
  } else if(mode == "current" || mode == "qcurrent" || mode == "ionic") {
    z = std::abs(track.GetDynamicParticle()->GetCharge()/CLHEP::eplus);
  } else {
    z = static_cast<G4double>(zProjectile);
  }

  G4int zInt = static_cast<G4int>(std::lround(z));
  zInt = std::max(1, zInt);

  if(!G4CCMottCorrection::IsAvailable(zInt) && G4CCMottCorrection::IsAvailable(zProjectile)) {
    zInt = zProjectile;
  }
  return zInt;
}

G4double G4CCBlendedDeltaChargeSquare(G4int Zp, G4int Qe, G4double deltaEnergy,
                                      G4double qEffective2, G4double globalWeight)
{
  const G4double w = std::max(0.0, std::min(1.0, globalWeight));
  const G4double q2Eff = std::max(0.0, qEffective2);
  if(w == 0.0) { return q2Eff; }
  const G4double z2Delta = G4CCScreening::DeltaChargeSquare(Zp, Qe, deltaEnergy);
  return std::max(0.0, w*z2Delta + (1.0 - w)*q2Eff);
}

std::atomic<int> gMottIntegralDiagnosticPrints{0};
std::atomic<int> gMottSampleDiagnosticPrints{0};

G4double G4CCDeltaIntegralMode(G4int Zp,
                               G4int zMott,
                               G4int Qe,
                               G4double tMin,
                               G4double tMax,
                               G4double beta2,
                               G4double totalEnergy2,
                               G4double spin,
                               G4double qEffective2,
                               G4double globalWeight,
                               G4bool useMott)
{
  if(tMax <= tMin || tMin <= 0.0 || beta2 <= 0.0 || Zp <= 0) { return 0.0; }

  const G4double beta = std::sqrt(beta2);
  constexpr G4int nBins = 160;
  const G4double u0 = std::log(tMin);
  const G4double u1 = std::log(tMax);
  const G4double du = (u1 - u0)/static_cast<G4double>(nBins);
  G4double integral = 0.0;

  for(G4int i = 0; i < nBins; ++i) {
    const G4double ua = u0 + du*static_cast<G4double>(i);
    const G4double ub = ua + du;
    const G4double um = 0.5*(ua + ub);
    const G4double ta = std::exp(ua);
    const G4double tm = std::exp(um);
    const G4double tb = std::exp(ub);
    const auto f = [&](G4double t) {
      const G4double z2 = G4CCBlendedDeltaChargeSquare(Zp, Qe, t, qEffective2, globalWeight);
      const G4double shape = useMott
        ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
        : G4CCMottCorrection::G4BornShape(beta2, t, tMax, totalEnergy2, spin);
      return z2*shape/t;
    };
    integral += (du/6.0)*(f(ta) + 4.0*f(tm) + f(tb));
  }
  return std::max(0.0, integral);
}
}

G4double G4CCIonDeltaElectronProcess::FormFactorThresholdEnergy(
  const G4ParticleDefinition* particle, G4double kineticEnergy)
{
  if(particle == nullptr || kineticEnergy <= 0.0) { return DBL_MAX; }
  const G4double eu = G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  if(G4CCDeltaFormFactorWeight(eu) <= 0.0) { return DBL_MAX; }
  const G4double ff = G4CCProjectileFormFactorParameter(particle, G4CCIonChargeState::GetZp(particle));
  return ff > 0.0 ? 1.e-6/ff : DBL_MAX;
}

G4double G4CCIonDeltaElectronProcess::FormFactorAcceptance(
  const G4ParticleDefinition* particle, G4double kineticEnergy, G4double deltaEnergy)
{
  if(particle == nullptr || kineticEnergy <= 0.0 || deltaEnergy <= 0.0) { return 1.0; }
  const G4double eu = G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  const G4double weight = G4CCDeltaFormFactorWeight(eu);
  if(weight <= 0.0) { return 1.0; }
  const G4double mass = particle->GetPDGMass();
  if(mass <= 0.0) { return 1.0; }
  const G4double totalEnergy = kineticEnergy + mass;
  const G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/(totalEnergy*totalEnergy);
  const G4double accepted = G4CCProjectileFormFactorAcceptance(deltaEnergy, beta2,
    MaxSecondaryEnergy(particle, kineticEnergy), totalEnergy*totalEnergy,
    particle->GetPDGSpin(), mass,
    G4CCProjectileFormFactorParameter(particle, G4CCIonChargeState::GetZp(particle)),
    G4CCProjectileMagMoment2(particle, mass));
  return G4CCClamp01(1.0 + weight*(accepted - 1.0));
}

G4CCIonDeltaElectronProcess::G4CCIonDeltaElectronProcess(const G4String& name)
  : G4VDiscreteProcess(name, fElectromagnetic),
    fAngularDistribution(new G4DeltaAngle())
{
  SetProcessSubType(92);
}

G4CCIonDeltaElectronProcess::~G4CCIonDeltaElectronProcess()
{
  delete fAngularDistribution;
}

G4bool G4CCIonDeltaElectronProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return G4CCIonChargeState::IsIonDefinition(&particle);
}

G4double G4CCIonDeltaElectronProcess::MaxSecondaryEnergy(
  const G4ParticleDefinition* particle, G4double kineticEnergy)
{
  if(particle == nullptr) { return 0.0; }
  const G4double mass = particle->GetPDGMass();
  if(mass <= 0.0 || kineticEnergy <= 0.0) { return 0.0; }
  const G4double tau = kineticEnergy/mass;
  const G4double ratio = electron_mass_c2/mass;
  const G4double denom = 1.0 + 2.0*(tau + 1.0)*ratio + ratio*ratio;
  if(denom <= 0.0) { return 0.0; }
  return 2.0*electron_mass_c2*tau*(tau + 2.0)/denom;
}

G4double G4CCIonDeltaElectronProcess::Shape(G4double deltaEnergy,
                                            G4double beta2,
                                            G4double tMax,
                                            G4double totalEnergy2,
                                            G4double spin)
{
  if(deltaEnergy <= 0.0 || tMax <= 0.0) { return 0.0; }
  G4double value = 1.0 - beta2*deltaEnergy/tMax;
  if(spin == 0.5 && totalEnergy2 > 0.0) {
    value += 0.5*deltaEnergy*deltaEnergy/totalEnergy2;
  }
  return std::max(0.0, value);
}

namespace
{
G4double G4CCFastBornShape(G4double deltaEnergy,
                              G4double beta2,
                              G4double tMax,
                              G4double totalEnergy2,
                              G4double spin)
{
  if(deltaEnergy <= 0.0 || tMax <= 0.0) { return 0.0; }
  G4double value = 1.0 - beta2*deltaEnergy/tMax;
  if(spin == 0.5 && totalEnergy2 > 0.0) {
    value += 0.5*deltaEnergy*deltaEnergy/totalEnergy2;
  }
  return std::max(0.0, value);
}

G4double G4CCAnalyticBornIntegralConstantZ2(G4double z2,
                                            G4double tMin,
                                            G4double tMax,
                                            G4double beta2,
                                            G4double totalEnergy2,
                                            G4double spin)
{
  if(z2 <= 0.0 || tMax <= tMin || tMin <= 0.0) { return 0.0; }

  G4double integral = 1.0/tMin - 1.0/tMax
    - (beta2/tMax)*std::log(tMax/tMin);
  if(spin == 0.5 && totalEnergy2 > 0.0) {
    integral += 0.5*(tMax - tMin)/totalEnergy2;
  }
  return std::max(0.0, z2*integral);
}

G4double G4CCGaussLogIntegral16(G4int Zp,
                                G4int zMott,
                                G4int Qe,
                                G4double tMin,
                                G4double tMax,
                                G4double beta2,
                                G4double totalEnergy2,
                                G4double spin,
                                G4double qEffective2,
                                G4double globalWeight,
                                G4bool useMott)
{
  static constexpr G4double x[16] = {
    -0.98940093499164993260, -0.94457502307323257608,
    -0.86563120238783174388, -0.75540440835500303390,
    -0.61787624440264374845, -0.45801677765722738634,
    -0.28160355077925891323, -0.09501250983763744019,
     0.09501250983763744019,  0.28160355077925891323,
     0.45801677765722738634,  0.61787624440264374845,
     0.75540440835500303390,  0.86563120238783174388,
     0.94457502307323257608,  0.98940093499164993260
  };
  static constexpr G4double w[16] = {
    0.02715245941175409485, 0.06225352393864789286,
    0.09515851168249278481, 0.12462897125553387205,
    0.14959598881657673208, 0.16915651939500253819,
    0.18260341504492358887, 0.18945061045506849629,
    0.18945061045506849629, 0.18260341504492358887,
    0.16915651939500253819, 0.14959598881657673208,
    0.12462897125553387205, 0.09515851168249278481,
    0.06225352393864789286, 0.02715245941175409485
  };

  const G4double beta = std::sqrt(beta2);
  const G4double u0 = std::log(tMin);
  const G4double u1 = std::log(tMax);
  const G4double mid = 0.5*(u0 + u1);
  const G4double half = 0.5*(u1 - u0);

  G4double integral = 0.0;
  for(G4int i = 0; i < 16; ++i) {
    const G4double t = std::exp(mid + half*x[i]);
    const G4double z2 = G4CCBlendedDeltaChargeSquare(Zp, Qe, t, qEffective2, globalWeight);
    const G4double shape = useMott
      ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
      : G4CCFastBornShape(t, beta2, tMax, totalEnergy2, spin);
    integral += w[i]*z2*shape/t;
  }
  return std::max(0.0, half*integral);
}
}

G4double G4CCIonDeltaElectronProcess::NumericalIntegral(G4int Zp,
                                                        G4int zMott,
                                                        G4int Qe,
                                                        G4double tMin,
                                                        G4double tMax,
                                                        G4double beta2,
                                                        G4double totalEnergy2,
                                                        G4double spin,
                                                        G4double qEffective2,
                                                        G4double globalWeight,
                                                        G4bool useMott)
{
  if(tMax <= tMin || tMin <= 0.0 || beta2 <= 0.0 || Zp <= 0) { return 0.0; }

  const G4double beta = std::sqrt(beta2);
  if(G4CCFastMottIntegralEnabled()) {
    const G4double wGlobal = std::max(0.0, std::min(1.0, globalWeight));

    if(Qe == 0 || wGlobal <= 0.0) {
      const G4double z2Const = (Qe == 0)
        ? G4CCBlendedDeltaChargeSquare(Zp, Qe, tMin, qEffective2, wGlobal)
        : std::max(0.0, qEffective2);
      if(useMott) {
        return std::max(0.0, z2Const
          * G4CCMottCorrection::IntegratedRMottOverT2(zMott, beta, tMin, tMax));
      }
      return G4CCAnalyticBornIntegralConstantZ2(z2Const, tMin, tMax, beta2,
                                                totalEnergy2, spin);
    }

    return G4CCGaussLogIntegral16(Zp, zMott, Qe, tMin, tMax, beta2, totalEnergy2,
                                  spin, qEffective2, wGlobal, useMott);
  }

  constexpr G4int nBins = 96;
  const G4double u0 = std::log(tMin);
  const G4double u1 = std::log(tMax);
  const G4double du = (u1 - u0)/static_cast<G4double>(nBins);
  G4double integral = 0.0;

  for(G4int i = 0; i < nBins; ++i) {
    const G4double ua = u0 + du*static_cast<G4double>(i);
    const G4double ub = ua + du;
    const G4double um = 0.5*(ua + ub);
    const G4double ta = std::exp(ua);
    const G4double tm = std::exp(um);
    const G4double tb = std::exp(ub);
    const auto f = [&](G4double t) {
      const G4double z2 = G4CCBlendedDeltaChargeSquare(Zp, Qe, t, qEffective2, globalWeight);
      const G4double shape = useMott
        ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
        : Shape(t, beta2, tMax, totalEnergy2, spin);
      return z2*shape/t;
    };
    integral += (du/6.0)*(f(ta) + 4.0*f(tm) + f(tb));
  }
  return std::max(0.0, integral);
}

G4double G4CCIonDeltaElectronProcess::ElectronProductionCut(const G4MaterialCutsCouple* couple)
{
  if(couple == nullptr) { return DBL_MAX; }
  const auto* cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
  if(cuts == nullptr) { return DBL_MAX; }
  const G4int idx = couple->GetIndex();
  if(idx < 0 || static_cast<std::size_t>(idx) >= cuts->size()) { return DBL_MAX; }
  return (*cuts)[static_cast<std::size_t>(idx)];
}

G4double G4CCIonDeltaElectronProcess::CrossSectionPerVolume(const G4Track& track) const
{
  const auto* particle = track.GetParticleDefinition();
  const auto* material = track.GetMaterial();
  if(!G4CCIonChargeState::IonisationModelApplies(particle, material)) { return 0.0; }

  const G4double kineticEnergy = track.GetKineticEnergy();
  const G4double tMax = MaxSecondaryEnergy(particle, kineticEnergy);
  const G4double tCut = ElectronProductionCut(track.GetMaterialCutsCouple());
  if(tMax <= tCut || tCut <= 0.0) { return 0.0; }

  const G4double mass = particle->GetPDGMass();
  const G4double totalEnergy = kineticEnergy + mass;
  const G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/(totalEnergy*totalEnergy);
  if(beta2 <= 0.0) { return 0.0; }

  const G4int Zp = G4CCIonChargeState::GetZp(particle);
  const G4int zMott = G4CCChooseMottZ(track, Zp);
  const G4int Qe = G4CCIonChargeState::GetQe(track);
  const G4double spin = particle->GetPDGSpin();
  const G4double eMeVPerU = G4CCIonChargeState::KineticEnergyPerNucleon(track)/CLHEP::MeV;
  if(!G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(eMeVPerU)) {
    return 0.0;
  }
  const G4double globalWeight = G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVPerU);
  const G4double qEffective2 = G4CCEffectiveChargeSquareForDelta(track);
  const G4bool useMott = G4CCUseMottDeltaElectronsAtEnergy(eMeVPerU)
    && G4CCMottCorrection::IsAvailable(zMott);
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) {
    G4CCDiracPhysics::CheckCut(tCut);
    const G4double native = globalWeight < 1. ?
      NumericalIntegral(Zp,zMott,Qe,tCut,tMax,beta2,totalEnergy*totalEnergy,spin,
                        qEffective2,0.,useMott)*twopi_mc2_rcl2/beta2 : 0.;
    const G4double tab = globalWeight > 0. ?
      G4CCDiracPhysics::Moments(particle,Qe,kineticEnergy,tCut)[0]*cm2 : 0.;
    return material->GetElectronDensity()*(globalWeight*tab+(1.-globalWeight)*native);
  }
  const G4double integral = NumericalIntegral(Zp, zMott, Qe, tCut, tMax, beta2,
                                              totalEnergy*totalEnergy, spin,
                                              qEffective2, globalWeight, useMott);

  if(G4CCMottDiagnosticsEnabled()) {
    const G4int n = gMottIntegralDiagnosticPrints.fetch_add(1);
    if(n < 20) {
      const G4double bornIntegral = G4CCDeltaIntegralMode(Zp, zMott, Qe, tCut, tMax, beta2,
        totalEnergy*totalEnergy, spin, qEffective2, globalWeight, false);
      const G4double mottIntegral = G4CCDeltaIntegralMode(Zp, zMott, Qe, tCut, tMax, beta2,
        totalEnergy*totalEnergy, spin, qEffective2, globalWeight, true);
      const G4double beta = std::sqrt(beta2);
      G4cout << "G4CC_MOTT_DIAGNOSTICS integral"
             << " Zp=" << Zp
             << " zMott=" << zMott
             << " zMottMode=" << G4CCMottZMode()
             << " Qe=" << Qe
             << " beta=" << beta
             << " tCut(MeV)=" << tCut/CLHEP::MeV
             << " tMax(MeV)=" << tMax/CLHEP::MeV
             << " useMott=" << useMott
             << " eMeVPerU=" << eMeVPerU
             << " mottMinMeVPerU=" << G4CCCurrentChargeConfig::MottMinimumEnergyMeVPerU()
             << " globalWeight=" << globalWeight
             << " qEffective2=" << qEffective2
             << " bornIntegral=" << bornIntegral
             << " mottIntegral=" << mottIntegral
             << " ratio=" << ((bornIntegral > 0.0) ? mottIntegral/bornIntegral : 0.0)
             << G4endl;
    }
  }

  if(integral <= 0.0) { return 0.0; }
  return material->GetElectronDensity()*twopi_mc2_rcl2*integral/beta2;
}

G4double G4CCIonDeltaElectronProcess::GetMeanFreePath(const G4Track& track,
                                                      G4double,
                                                      G4ForceCondition* condition)
{
  if(condition != nullptr) { *condition = NotForced; }
  const G4double sigma = CrossSectionPerVolume(track);
  if(sigma <= 0.0) { return DBL_MAX; }
  return 1.0/sigma;
}

G4double G4CCIonDeltaElectronProcess::SampleDeltaEnergy(const G4Track& track,
                                                        G4double tMin,
                                                        G4double tMax,
                                                        G4double beta2,
                                                        G4double totalEnergy2,
                                                        G4double spin,
                                                        G4double qEffective2,
                                                        G4double globalWeight,
                                                        G4bool useMott) const
{
  const G4int Zp = G4CCIonChargeState::GetZp(track.GetParticleDefinition());
  const G4int zMott = G4CCChooseMottZ(track, Zp);
  const G4int Qe = G4CCIonChargeState::GetQe(track);
  if(Zp <= 0 || tMax <= tMin || tMin <= 0.0) { return 0.0; }

  const G4double bareZ2 = static_cast<G4double>(Zp)*static_cast<G4double>(Zp);
  const G4double envelopeZ2 = std::max(bareZ2, std::max(1.0, qEffective2));
  const G4double beta = std::sqrt(beta2);

  G4double shapeMax = useMott
    ? G4CCMottCorrection::MaxRMott(zMott, beta, tMin, tMax)
    : 1.0;
  if(!useMott && spin == 0.5 && totalEnergy2 > 0.0) {
    shapeMax += 0.5*tMax*tMax/totalEnergy2;
  }
  shapeMax = std::max(shapeMax, 1.0);

  for(G4int tries = 0; tries < 10000; ++tries) {
    const G4double q = G4UniformRand();
    const G4double t = tMin*tMax/(tMin*(1.0 - q) + tMax*q);
    const G4double z2 = std::min(envelopeZ2,
      G4CCBlendedDeltaChargeSquare(Zp, Qe, t, qEffective2, globalWeight));
    const G4double shape = useMott
      ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
      : Shape(t, beta2, tMax, totalEnergy2, spin);
    const G4double accept = (z2/envelopeZ2)*shape/shapeMax;
    if(G4UniformRand() <= accept) { return t; }
  }

  return 0.0;
}

G4VParticleChange* G4CCIonDeltaElectronProcess::PostStepDoIt(const G4Track& track,
                                                             const G4Step&)
{
  ClearNumberOfInteractionLengthLeft();

  fParticleChange.InitializeForPostStep(track);

  const auto* particle = track.GetParticleDefinition();
  const auto* material = track.GetMaterial();
  if(!G4CCIonChargeState::IonisationModelApplies(particle, material)) {
    return &fParticleChange;
  }

  const G4double kineticEnergy = track.GetKineticEnergy();
  const G4double mass = particle->GetPDGMass();
  const G4double totalEnergy = kineticEnergy + mass;
  const G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/(totalEnergy*totalEnergy);
  const G4double tMax = MaxSecondaryEnergy(particle, kineticEnergy);
  const G4double tCut = ElectronProductionCut(track.GetMaterialCutsCouple());
  if(beta2 <= 0.0 || tMax <= tCut || tCut <= 0.0) { return &fParticleChange; }

  const G4double eMeVPerU = G4CCIonChargeState::KineticEnergyPerNucleon(track)/CLHEP::MeV;
  if(!G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(eMeVPerU)) {
    return &fParticleChange;
  }
  const G4double globalWeight = G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVPerU);
  const G4double qEffective2 = G4CCEffectiveChargeSquareForDelta(track);
  const G4int Zp = G4CCIonChargeState::GetZp(particle);
  const G4double formFactorWeight = G4CCDeltaFormFactorWeight(eMeVPerU);
  const G4double projectileFormFactor = G4CCProjectileFormFactorParameter(particle, Zp);
  const G4double projectileMagMoment2 = G4CCProjectileMagMoment2(particle, mass);
  const G4int zMottForGate = G4CCChooseMottZ(track, Zp);
  const G4bool useMott = G4CCUseMottDeltaElectronsAtEnergy(eMeVPerU)
    && G4CCMottCorrection::IsAvailable(zMottForGate);

  G4double deltaEnergy = 0.;
  G4bool sampledDirac = false;
  const G4bool diracMode=G4CCCurrentChargeConfig::Instance().IsDiracMode();
  if(diracMode) {
    G4CCDiracPhysics::CheckCut(tCut);
    const auto n=G4CCIonChargeState::GetQe(track);
    const G4double tab = globalWeight>0. ? globalWeight*
      G4CCDiracPhysics::Moments(particle,n,kineticEnergy,tCut)[0]*cm2 : 0.;
    const G4double native = globalWeight<1. ? (1.-globalWeight)*
      NumericalIntegral(Zp,zMottForGate,n,tCut,tMax,beta2,totalEnergy*totalEnergy,
                        particle->GetPDGSpin(),qEffective2,0.,useMott)*twopi_mc2_rcl2/beta2 : 0.;
    if(tab+native<=0.) return &fParticleChange;
    sampledDirac=G4UniformRand()*(tab+native)<tab;
    deltaEnergy=sampledDirac ?
      G4CCDiracPhysics::Sample(particle,n,kineticEnergy,tCut,G4UniformRand(),G4UniformRand()) :
      SampleDeltaEnergy(track,tCut,tMax,beta2,totalEnergy*totalEnergy,
                        particle->GetPDGSpin(),qEffective2,0.,useMott);
  } else {
    deltaEnergy=SampleDeltaEnergy(track,tCut,tMax,beta2,totalEnergy*totalEnergy,
                                 particle->GetPDGSpin(),qEffective2,globalWeight,useMott);
  }
  if(G4CCMottDiagnosticsEnabled() && deltaEnergy > 0.0) {
    const G4int n = gMottSampleDiagnosticPrints.fetch_add(1);
    if(n < 20) {
      const G4int Zp = G4CCIonChargeState::GetZp(particle);
      const G4int zMott = G4CCChooseMottZ(track, Zp);
      const G4double beta = std::sqrt(beta2);
      const G4double born = G4CCMottCorrection::G4BornShape(beta2, deltaEnergy, tMax,
        totalEnergy*totalEnergy, particle->GetPDGSpin());
      const G4double mott = G4CCMottCorrection::RMott(zMott, beta, deltaEnergy, tMax);
      G4cout << "G4CC_MOTT_DIAGNOSTICS sample"
             << " Zp=" << Zp
             << " zMott=" << zMott
             << " zMottMode=" << G4CCMottZMode()
             << " useMott=" << useMott
             << " eMeVPerU=" << eMeVPerU
             << " mottMinMeVPerU=" << G4CCCurrentChargeConfig::MottMinimumEnergyMeVPerU()
             << " globalWeight=" << globalWeight
             << " deltaE(MeV)=" << deltaEnergy/CLHEP::MeV
             << " bornShape=" << born
             << " mottShape=" << mott
             << " C_Mott=" << ((born > 0.0) ? mott/born : 0.0)
             << G4endl;
    }
  }
  if(deltaEnergy <= 0.0 || deltaEnergy >= kineticEnergy) { return &fParticleChange; }

  const G4double ffAccept = G4CCProjectileFormFactorAcceptance(deltaEnergy, beta2, tMax,
    totalEnergy*totalEnergy, particle->GetPDGSpin(), mass, projectileFormFactor,
    projectileMagMoment2);
  const G4double mixAccept = diracMode ?
    (sampledDirac ? 1. : 1.-formFactorWeight+formFactorWeight*ffAccept) :
    G4CCFallbackFormFactorAcceptance(Zp, G4CCIonChargeState::GetQe(track),
      deltaEnergy, qEffective2, globalWeight, formFactorWeight, ffAccept);
  if(G4UniformRand() > mixAccept) { return &fParticleChange; }

  const G4double deltaMomentum = std::sqrt(deltaEnergy*(deltaEnergy + 2.0*electron_mass_c2));
  const G4double primaryMomentum = totalEnergy*std::sqrt(std::max(0.0, beta2));
  if(deltaMomentum <= 0.0 || primaryMomentum <= 0.0) { return &fParticleChange; }

  G4ThreeVector deltaDirection;
  if(G4CCUseDeltaAngularGeneratorRuntime() && fAngularDistribution != nullptr) {
    const G4int zTarget = G4CCSelectRandomTargetZ(material);
    if(zTarget > 0) {
      deltaDirection = fAngularDistribution->SampleDirection(track.GetDynamicParticle(),
                                                             deltaEnergy, zTarget, material);
    }
  }

  if(deltaDirection.mag2() <= 0.0) {
    G4double cosTheta = deltaEnergy*(totalEnergy + electron_mass_c2)/(deltaMomentum*primaryMomentum);
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
    const G4double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta)*(1.0 + cosTheta)));
    const G4double phi = twopi*G4UniformRand();

    deltaDirection.set(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
    deltaDirection.rotateUz(track.GetMomentumDirection());
  }

  auto* electron = new G4DynamicParticle(G4Electron::Electron(), deltaDirection, deltaEnergy);
  auto* secondaryTrack = new G4Track(electron, track.GetGlobalTime(), track.GetPosition());
  secondaryTrack->SetTouchableHandle(track.GetTouchableHandle());
  fParticleChange.SetNumberOfSecondaries(1);
  fParticleChange.AddSecondary(secondaryTrack);

  const G4double newKineticEnergy = std::max(0.0, kineticEnergy - deltaEnergy);
  const G4ThreeVector initialMomentum = track.GetMomentumDirection()*primaryMomentum;
  const G4ThreeVector finalMomentum = initialMomentum - deltaDirection*deltaMomentum;
  if(finalMomentum.mag2() > 0.0) {
    fParticleChange.SetProposedMomentumDirection(finalMomentum.unit());
  }
  fParticleChange.SetProposedKineticEnergy(newKineticEnergy);
  fParticleChange.SetProposedCharge(track.GetDynamicParticle()->GetCharge());

  return &fParticleChange;
}

void G4CCIonDeltaElectronProcess::PreparePhysicsTable(const G4ParticleDefinition& particle)
{
  G4VDiscreteProcess::PreparePhysicsTable(particle);
  if(!G4CCCurrentChargeConfig::Instance().IsDiracMode()) return;
  G4CCDiracPhysics::Tables();
  G4CCGlobalChargeExchangeTable::Instance().EnsureLoaded();
  if(!G4CCGlobalChargeExchangeTable::Instance().IsLoaded()) {
    G4Exception("CCDiracPrepare","CCDiracGLOBAL",FatalException,"Required GLOBAL data did not load.");
  }
  G4CCAtimaEquilibriumDedxTable::Instance();
  auto* table=G4ProductionCutsTable::GetProductionCutsTable();
  const auto* cuts=table->GetEnergyCutsVector(1);
  if(!cuts || cuts->size()<table->GetTableSize()) {
    G4Exception("CCDiracPrepare","CCDiracCuts",FatalException,
                "Electron production cuts were not initialized for all material/region couples.");
    return;
  }
  for(std::size_t i=0;i<table->GetTableSize();++i) {
    const auto* couple=table->GetMaterialCutsCouple(i);
    const auto* mat=couple->GetMaterial();
    if(!G4CCDiracPhysics::MaterialInDomain(mat))continue;
    bool applicable=false;
    for(int z=29;z<=92;++z) {
      if(!G4CCDiracPhysics::Tables().HasProjectile(z))continue;
      G4CCDiracPhysics::RequireMaterialCoverage(z,mat);
      for(const auto* el:*mat->GetElementVector())
        G4CCAtimaEquilibriumDedxTable::Instance().RequirePair(z,el->GetZasInt());
      applicable=true;
    }
    if(applicable) G4CCDiracPhysics::CheckCut((*cuts)[i]);
  }
}
