#include "G4CCNonEMChargeResetProcess.hh"

#include "G4CCIonChargeState.hh"
#include "G4ForceCondition.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

G4CCNonEMChargeResetProcess::G4CCNonEMChargeResetProcess(const G4String& name)
  : G4VDiscreteProcess(name, fUserDefined)
{
  SetProcessSubType(93);
}

G4bool G4CCNonEMChargeResetProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return G4CCIonChargeState::IsIonDefinition(&particle);
}

G4double G4CCNonEMChargeResetProcess::PostStepGetPhysicalInteractionLength(
  const G4Track&, G4double, G4ForceCondition* condition)
{
  if(condition != nullptr) { *condition = Forced; }
  return DBL_MAX;
}

G4double G4CCNonEMChargeResetProcess::GetMeanFreePath(const G4Track&,
                                                      G4double,
                                                      G4ForceCondition* condition)
{
  if(condition != nullptr) { *condition = Forced; }
  return DBL_MAX;
}

G4bool G4CCNonEMChargeResetProcess::ShouldResetAfter(const G4VProcess* process)
{
  if(process == nullptr) { return false; }
  const G4ProcessType type = process->GetProcessType();
  if(type == fTransportation || type == fElectromagnetic || type == fParallel) { return false; }

  const G4String& name = process->GetProcessName();
  static const std::vector<G4String> excluded = {
    "Transportation", "CoupledTransportation", "StepLimiter", "UserSpecialCuts", "ionElastic",
    "CCNonEMChargeReset", "CCIonChargeExchange", "CCIonDeltaElectron", "CCIonIoni"
  };
  return std::find(excluded.begin(), excluded.end(), name) == excluded.end();
}

G4VParticleChange* G4CCNonEMChargeResetProcess::PostStepDoIt(const G4Track& track,
                                                             const G4Step& step)
{
  ClearNumberOfInteractionLengthLeft();

  fParticleChange.InitializeForPostStep(track);
  if(!G4CCIonChargeState::CurrentChargeModelApplies(track)) { return &fParticleChange; }

  const G4int Zp = G4CCIonChargeState::GetZp(track);
  if(Zp <= 0) { return &fParticleChange; }

  const auto* post = step.GetPostStepPoint();
  const auto* process = post == nullptr ? nullptr : post->GetProcessDefinedStep();
  if(!ShouldResetAfter(process)) { return &fParticleChange; }

  G4CCIonChargeState::ProposeBareIon(fParticleChange, Zp);
  return &fParticleChange;
}
