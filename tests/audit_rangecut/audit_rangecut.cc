#include "G4CCPhysicsOptions.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4Box.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4Positron.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
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
std::size_t checks=0, rows=0;
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
    auto* vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    silicon=G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    auto* w=new G4LogicalVolume(new G4Box("world",10*mm,10*mm,10*mm),vacuum,"world");
    auto* pv=new G4PVPlacement(nullptr,{},w,"world",nullptr,false,0);
    auto* s=new G4LogicalVolume(new G4Box("silicon",5*mm,5*mm,.15*mm),silicon,"silicon");
    new G4PVPlacement(nullptr,{},s,"silicon",w,false,0);
    return pv;
  }
};
class Physics final:public G4VModularPhysicsList {
 public:
  Physics() {
    SetVerboseLevel(0);SetDefaultCutValue(25*um);
    RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(options,0));
  }
  void SetCuts() override {SetCutsWithDefault();}
};
class Gun final:public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{1};
  void GeneratePrimaries(G4Event* event) override {
    gun.SetParticleDefinition(G4Electron::Electron());
    gun.SetParticleEnergy(50*MeV);
    gun.SetParticlePosition({0,0,-1*mm});
    gun.SetParticleMomentumDirection({0,0,1});
    gun.GeneratePrimaryVertex(event);
  }
};
}

int main(int argc,char** argv) {
  try {
    if(argc!=6)throw std::runtime_error("Usage: audit_rangecut delta stopping global atima-low output.csv");
    options.tables.delta=argv[1];options.tables.stopping=argv[2];
    options.tables.global=argv[3];options.tables.atimaLow=argv[4];
    options.scattering.mode=G4CCScatteringMode::Opt4MultipleScattering;
    Exceptions exceptions;
    G4StateManager::GetStateManager()->SetExceptionHandler(&exceptions);
    auto run=std::make_unique<G4RunManager>();
    run->SetUserInitialization(new Detector);
    run->SetUserInitialization(new Physics);
    run->SetUserAction(new Gun);
    run->Initialize();run->BeamOn(0);
    auto* ion=G4IonTable::GetIonTable()->GetIon(54,124,0);
    std::ofstream out(argv[5]);
    Require(bool(out),"Cannot create output");
    out<<std::setprecision(17)<<"particle,kinetic_MeV,sequence_index,range_cut_mm,gamma_cut_keV,electron_cut_keV,positron_cut_keV,range_api_MeV_per_mm,explicit_cut_MeV_per_mm\n";
    G4EmCalculator reused;
    reused.SetVerbose(0);
    const std::array<double,8> sequence{{1*um,1*mm,1*um,10*mm,25*um,1*mm,100*um,1*um}};
    const std::array<double,4> ionEu{{35.,400.,1750.,2500.}};
    double cutFirst=0.,cutSecond=0.,lossFirst=0.,lossSecond=0.;
    const auto test=[&](const G4ParticleDefinition* particle,double kinetic) {
      for(std::size_t index=0;index<sequence.size();++index) {
        const double range=sequence[index];
        G4EmCalculator fresh;
        fresh.SetVerbose(0);
        const double gammaCut=fresh.ComputeEnergyCutFromRangeCut(range,G4Gamma::Gamma(),silicon);
        const double electronCut=fresh.ComputeEnergyCutFromRangeCut(range,G4Electron::Electron(),silicon);
        const double positronCut=fresh.ComputeEnergyCutFromRangeCut(range,G4Positron::Positron(),silicon);
        Require(gammaCut>0.&&electronCut>0.&&positronCut>0.,"Invalid range-to-energy conversion");
        double expected=0.;
        if(particle==ion) {
          expected=fresh.ComputeDEDX(kinetic,particle,"CCIonIoni",silicon,electronCut);
        } else {
          expected=fresh.ComputeDEDX(kinetic,particle,"eIoni",silicon,electronCut)
            +fresh.ComputeDEDX(kinetic,particle,"eBrem",silicon,gammaCut);
        }
        const double actual=reused.ComputeDEDXForCutInRange(kinetic,particle,silicon,range);
        Require(expected>0.,"Nonpositive explicit-cut reference");
        out<<particle->GetParticleName()<<','<<kinetic/MeV<<','<<index<<','<<range/mm<<','
           <<gammaCut/keV<<','<<electronCut/keV<<','<<positronCut/keV<<','
           <<actual/(MeV/mm)<<','<<expected/(MeV/mm)<<'\n';
        out.flush();
        ++rows;
        Equal(actual,expected,"Range API disagrees with fresh explicit-cut reference at index "+std::to_string(index));
        if(particle==ion&&std::abs(kinetic/ion->GetPDGMass()*931.494-400.)<1.e-8) {
          if(index==0){cutFirst=electronCut;lossFirst=actual;}
          if(index==1){cutSecond=electronCut;lossSecond=actual;}
          if(index==2||index==7)Equal(actual,lossFirst,"Returning to 1 um must restore original stopping");
        }
      }
    };
    for(double eu:ionEu)test(ion,eu*MeV*ion->GetPDGMass()/(931.494*MeV));
    test(G4Electron::Electron(),50*MeV);
    test(G4Positron::Positron(),50*MeV);
    Require(cutSecond>100.*cutFirst,"Test did not exercise substantially different cuts");
    Require(lossSecond>1.2*lossFirst,"Xe restricted stopping did not respond to the larger cut");
    std::cout<<std::setprecision(17)<<"AUDIT_RANGECUT_PASS checks="<<checks<<" rows="<<rows
      <<" cut_1um_keV="<<cutFirst/keV<<" cut_1mm_keV="<<cutSecond/keV
      <<" Xe400_1um_MeV_per_mm="<<lossFirst/(MeV/mm)
      <<" Xe400_1mm_MeV_per_mm="<<lossSecond/(MeV/mm)<<'\n';
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"AUDIT_RANGECUT_FAIL "<<e.what()<<'\n';return 1;
  }
}
