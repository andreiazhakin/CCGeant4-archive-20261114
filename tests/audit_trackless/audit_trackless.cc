#include "G4CCPhysicsOptions.hh"
#include "G4CCLindhardSorensenIonModel.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4Box.hh"
#include "G4DynamicParticle.hh"
#include "G4EmCalculator.hh"
#include "G4Event.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VExceptionHandler.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
G4CCPhysicsOptions options;
G4Material* silicon=nullptr;
G4Material* vacuum=nullptr;
std::size_t checks=0;

class Exceptions final:public G4VExceptionHandler {
  G4bool Notify(const char* origin,const char* code,G4ExceptionSeverity severity,
                const char* description) override {
    std::cerr<<origin<<' '<<code<<' '<<description<<'\n';
    if(severity==FatalException||severity==FatalErrorInArgument)
      throw std::runtime_error(std::string(code)+": "+description);
    return false;
  }
};

void Require(bool ok,const std::string& message) {
  ++checks;
  if(!ok)throw std::runtime_error(message);
}
void Equal(double a,double b,const std::string& message) {
  Require(std::isfinite(a)&&std::isfinite(b)&&
    std::abs(a-b)<=2.e-12*std::max({1.e-30,std::abs(a),std::abs(b)}),message);
}

class Detector final:public G4VUserDetectorConstruction {
  G4VPhysicalVolume* Construct() override {
    vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    silicon=G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    auto* w=new G4LogicalVolume(new G4Box("world",10*mm,10*mm,10*mm),vacuum,"world");
    auto* pv=new G4PVPlacement(nullptr,{},w,"world",nullptr,false,0);
    auto* s=new G4LogicalVolume(new G4Box("silicon",5*mm,5*mm,0.15*mm),silicon,"silicon");
    new G4PVPlacement(nullptr,{},s,"silicon",w,false,0);
    return pv;
  }
};
class Physics final:public G4VModularPhysicsList {
 public:
  Physics() {
    SetVerboseLevel(0);
    SetDefaultCutValue(25*um);
    RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(options,0));
  }
  void SetCuts() override {SetCutsWithDefault();}
};
class Gun final:public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{1};
  void GeneratePrimaries(G4Event* event) override {
    auto* ion=G4IonTable::GetIonTable()->GetIon(54,124,0);
    gun.SetParticleDefinition(ion);
    gun.SetParticleEnergy(400*MeV*ion->GetPDGMass()/(931.494*MeV));
    gun.SetParticlePosition({0,0,-1*mm});
    gun.SetParticleMomentumDirection({0,0,1});
    gun.GeneratePrimaryVertex(event);
  }
};

void Prime(G4CCLindhardSorensenIonModel& model,const G4ParticleDefinition* ion,
           const G4MaterialCutsCouple* couple,int electrons,double kinetic) {
  G4Step step;
  auto* dp=new G4DynamicParticle(ion,{0,0,1},kinetic);
  dp->SetCharge((54-electrons)*eplus);
  G4Track track(dp,0,{});
  track.SetTrackID(100+electrons);
  track.SetStep(&step);
  step.SetTrack(&track);
  step.GetPreStepPoint()->SetMaterial(const_cast<G4Material*>(couple->GetMaterial()));
  step.GetPreStepPoint()->SetMaterialCutsCouple(couple);
  model.ChargeSquareRatio(track);
  model.SetCurrentCouple(couple);
}
}

int main(int argc,char** argv) {
  try {
    if(argc!=6)throw std::runtime_error("Usage: audit_trackless delta stopping global atima-low output.csv");
    options.tables.delta=argv[1];options.tables.stopping=argv[2];
    options.tables.global=argv[3];options.tables.atimaLow=argv[4];
    options.scattering.mode=G4CCScatteringMode::Opt4MultipleScattering;
    Exceptions exceptions;
    G4StateManager::GetStateManager()->SetExceptionHandler(&exceptions);
    auto run=std::make_unique<G4RunManager>();
    run->SetUserInitialization(new Detector);
    run->SetUserInitialization(new Physics);
    run->SetUserAction(new Gun);
    run->Initialize();
    run->BeamOn(0);
    auto* ion=G4IonTable::GetIonTable()->GetIon(54,124,0);
    auto* pm=G4GenericIon::GenericIon()->GetProcessManager();
    G4CCLindhardSorensenIonModel* model=nullptr;
    std::string process;
    for(G4int i=0;i<pm->GetProcessListLength();++i) {
      auto* loss=dynamic_cast<G4VEnergyLossProcess*>((*pm->GetProcessList())[i]);
      if(loss) {
        auto* candidate=dynamic_cast<G4CCLindhardSorensenIonModel*>(loss->GetModelByIndex(0));
        if(candidate){model=candidate;process=loss->GetProcessName();break;}
      }
    }
    Require(model!=nullptr,"No initialized CC Lindhard-Sorensen model found");
    const G4MaterialCutsCouple* siCouple=nullptr;
    const G4MaterialCutsCouple* vacuumCouple=nullptr;
    const auto* cuts=G4ProductionCutsTable::GetProductionCutsTable();
    for(std::size_t i=0;i<cuts->GetTableSize();++i) {
      const auto* c=cuts->GetMaterialCutsCouple(i);
      if(c->GetMaterial()==silicon)siCouple=c;
      if(c->GetMaterial()==vacuum)vacuumCouple=c;
    }
    Require(siCouple&&vacuumCouple,"Missing initialized material couples");
    const double materialCut=(*cuts->GetEnergyCutsVector(1))[siCouple->GetIndex()];
    Require(std::abs(materialCut-keV)>1.e-6*keV,"Production cut must differ from caller cut");
    std::ofstream out(argv[5]);
    Require(bool(out),"Cannot create output");
    out<<std::setprecision(17)<<"Eu_MeV_per_u,cut_keV,previous_N,previous_material,expected_bare_MeV_per_mm,calculator_MeV_per_mm,direct_MeV_per_mm\n";
    G4EmCalculator calculator;
    calculator.SetVerbose(0);
    const std::array<double,7> energies{{1.,20.,35.,400.,1000.,1750.,2500.}};
    const std::array<double,4> cutValues{{.99,1.,10.,100.}};
    const std::array<int,4> histories{{0,1,15,27}};
    const std::array<const G4MaterialCutsCouple*,2> couples{{siCouple,vacuumCouple}};
    for(double eu:energies)for(double ck:cutValues) {
      const double kinetic=eu*MeV*ion->GetPDGMass()/(931.494*MeV);
      const double cut=ck*keV;
      double expected=0.;
      const bool explicitResult=model->ComputeDiracRestrictedDEDX(silicon,ion,kinetic,cut,54.,expected);
      if(!explicitResult)expected=calculator.ComputeDEDX(kinetic,ion,process,silicon,cut);
      Require(expected>0.,"Nonpositive reference stopping");
      for(const auto* previous:couples)for(int n:histories) {
        Prime(*model,ion,previous,n,400*MeV*ion->GetPDGMass()/(931.494*MeV));
        const double actual=calculator.ComputeDEDX(kinetic,ion,process,silicon,cut);
        Equal(actual,expected,"Real calculator depends on previous charge/material/cut");
        Prime(*model,ion,previous,n,400*MeV*ion->GetPDGMass()/(931.494*MeV));
        model->GetChargeSquareRatio(ion,silicon,kinetic);
        const double length=nm;
        double eloss=length*model->ComputeDEDXPerVolume(silicon,ion,kinetic,cut);
        model->CorrectionsAlongStep(silicon,ion,kinetic,cut,length,eloss);
        Equal(eloss/length,expected,"Direct calculator sequence depends on previous track");
        if(eu>30.)Equal(model->GetParticleCharge(ion,silicon,kinetic)/eplus,54.,
                       "Trackless charge getter retained transport cache");
        out<<eu<<','<<ck<<','<<n<<','<<previous->GetMaterial()->GetName()<<','
           <<expected/(MeV/mm)<<','<<actual/(MeV/mm)<<','<<eloss/length/(MeV/mm)<<'\n';
      }
    }
    const double k400=400*MeV*ion->GetPDGMass()/(931.494*MeV);
    double n0=0.,n15=0.,n27=0.;
    Require(model->ComputeDiracRestrictedDEDX(silicon,ion,k400,keV,54.,n0),"No explicit N=0 stopping");
    Require(model->ComputeDiracRestrictedDEDX(silicon,ion,k400,keV,39.,n15),"No explicit N=15 stopping");
    Require(model->ComputeDiracRestrictedDEDX(silicon,ion,k400,keV,27.,n27),"No explicit N=27 stopping");
    Require(n0>n15&&n15>n27,"Explicit-charge API lost state dependence");
    std::cout<<std::setprecision(17)<<"AUDIT_TRACKLESS_PASS checks="<<checks
      <<" Si_production_cut_keV="<<materialCut/keV<<" N0="<<n0/(MeV/mm)
      <<" N15="<<n15/(MeV/mm)<<" N27="<<n27/(MeV/mm)<<'\n';
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"AUDIT_TRACKLESS_FAIL "<<e.what()<<'\n';
    return 1;
  }
}
