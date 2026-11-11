#include "G4CCPhysicsOptions.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCDiracPhysics.hh"
#include "G4CCDiracTable.hh"
#include "G4CCIonChargeState.hh"
#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCIonDeltaElectronProcess.hh"
#include "G4CCLindhardSorensenIonModel.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4Box.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4EmParameters.hh"
#include "G4Event.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MTRunManager.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VProcess.hh"
#include "G4VExceptionHandler.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double kAtomicMassMeV = 931.494;
struct Settings {
  int Z=54, A=124, N=-1, events=400, threads=1;
  long seed=6110101;
  double eu=400, thickness=100, cut=25, cutKeV=0, secondaryEu=0;
  std::string delta, stopping, global, atimaLow, ss, out, material="G4_Si";
  std::string scattering="msc";
  bool kernel=true;
  bool legacy=false;
  bool trace=false,secondary=false;
  std::array<double,4> edges{{30,40,1500,2000}};
} cfg;
std::mutex resultsMutex;
std::ofstream output;
std::ofstream stepOutput;
G4Material* slabMaterial=nullptr;
std::size_t completed=0;
double maxPrimaryResidual=0;
const auto start=std::chrono::steady_clock::now();

class TestExceptionHandler final:public G4VExceptionHandler {
  G4bool Notify(const char* origin,const char* code,G4ExceptionSeverity severity,
                const char* description) override {
    std::cerr<<"G4_EXCEPTION "<<origin<<' '<<code<<' '<<description<<'\n';
    if(severity==FatalException||severity==FatalErrorInArgument)
      throw std::runtime_error(std::string(code)+": "+description);
    return false;
  }
};

void Parse(int argc,char** argv) {
  if((argc-1)%2) throw std::runtime_error("Expected --key value pairs");
  for(int i=1;i<argc;i+=2) {
    const std::string k=argv[i],v=argv[i+1];
    if(k=="--delta")cfg.delta=v;
    else if(k=="--stopping")cfg.stopping=v;
    else if(k=="--global")cfg.global=v;
    else if(k=="--atima-low")cfg.atimaLow=v;
    else if(k=="--ss")cfg.ss=v;
    else if(k=="--scatter")cfg.scattering=v;
    else if(k=="--kernel")cfg.kernel=std::stoi(v)!=0;
    else if(k=="--legacy")cfg.legacy=std::stoi(v)!=0;
    else if(k=="--trace")cfg.trace=std::stoi(v)!=0;
    else if(k=="--secondary")cfg.secondary=std::stoi(v)!=0;
    else if(k=="--secondary-energy")cfg.secondaryEu=std::stod(v);
    else if(k=="--material")cfg.material=v;
    else if(k=="--out")cfg.out=v;
    else if(k=="--Z")cfg.Z=std::stoi(v);
    else if(k=="--A")cfg.A=std::stoi(v);
    else if(k=="--N")cfg.N=std::stoi(v);
    else if(k=="--events")cfg.events=std::stoi(v);
    else if(k=="--threads")cfg.threads=std::stoi(v);
    else if(k=="--seed")cfg.seed=std::stol(v);
    else if(k=="--energy")cfg.eu=std::stod(v);
    else if(k=="--thickness")cfg.thickness=std::stod(v);
    else if(k=="--cut")cfg.cut=std::stod(v);
    else if(k=="--cut-keV")cfg.cutKeV=std::stod(v);
    else if(k=="--low0")cfg.edges[0]=std::stod(v);
    else if(k=="--low1")cfg.edges[1]=std::stod(v);
    else if(k=="--high0")cfg.edges[2]=std::stod(v);
    else if(k=="--high1")cfg.edges[3]=std::stod(v);
    else throw std::runtime_error("Unknown argument: "+k);
  }
  if(cfg.out.empty()||cfg.events<0||cfg.threads<1||cfg.Z<1||cfg.A<cfg.Z||
     cfg.N < -1 || cfg.N>cfg.Z||cfg.eu<=0||cfg.thickness<=0||cfg.cut<=0||cfg.cutKeV<0||cfg.secondaryEu<0)
    throw std::runtime_error("Invalid test configuration");
  if(cfg.scattering!="msc"&&cfg.scattering!="native"&&cfg.scattering!="table")
    throw std::runtime_error("--scatter must be msc, native, or table");
}

class Detector final:public G4VUserDetectorConstruction {
  G4VPhysicalVolume* Construct() override {
    auto* nist=G4NistManager::Instance();
    auto* vacuum=nist->FindOrBuildMaterial("G4_Galactic");
    auto* worldLV=new G4LogicalVolume(new G4Box("World",10*mm,10*mm,10*mm+cfg.thickness*um),vacuum,"World");
    auto* world=new G4PVPlacement(nullptr,{},worldLV,"World",nullptr,false,0,true);
    slabMaterial=nist->FindOrBuildMaterial(cfg.material);
    auto* slabLV=new G4LogicalVolume(new G4Box("Slab",5*mm,5*mm,cfg.thickness*.5*um),slabMaterial,"Slab");
    new G4PVPlacement(nullptr,{},slabLV,"Slab",worldLV,false,0,true);
    return world;
  }
};

class Physics final:public G4VModularPhysicsList {
public:
  Physics() {
    SetVerboseLevel(0);
    SetDefaultCutValue(cfg.cut*um);
    G4CCPhysicsOptions opts;
    opts.tables.delta=cfg.delta;
    opts.tables.stopping=cfg.stopping;
    opts.tables.global=cfg.global;
    opts.tables.atimaLow=cfg.atimaLow;
    opts.transition={cfg.edges[0]*MeV,cfg.edges[1]*MeV,cfg.edges[2]*MeV,cfg.edges[3]*MeV};
    opts.scattering.mode=cfg.scattering=="msc"?G4CCScatteringMode::Opt4MultipleScattering:
      (cfg.scattering=="native"?G4CCScatteringMode::NativeSingleScattering:G4CCScatteringMode::TabulatedSingleScattering);
    opts.scattering.tables=cfg.ss;
    if(cfg.legacy) {
      using Mode=G4EmStandardPhysics_option4_CurrentCharge::IonScatteringMode;
      const auto mode=cfg.scattering=="msc"?Mode::Opt4MultipleScattering:
        (cfg.scattering=="native"?Mode::NativeSingleScattering:Mode::TabulatedSingleScattering);
      RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(0,mode));
    } else RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(opts,0));
  }
  void SetCuts() override {
    SetCutsWithDefault();
    if(cfg.cutKeV>0) {
      G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyCutVector(
        std::vector<G4double>(64,cfg.cutKeV*keV),1);
    }
  }
};

struct EventData {
  int event=0,N=0,deltaCount=0,chargeExchange=0,scatter=0,primarySteps=0;
  double initial=0,initialQ=0,edepSlab=0,edepAll=0,deltaEnergy=0;
  double finalK=0,finalQ=0,finalAngle=0,primaryResidual=0,maxStepResidual=0;
  double escapedK=0,maxDeltaRatio=0;
  double injectedK=0,secondaryInitialQ=0,secondaryFinalQ=0,secondaryFinalK=0;
  int syntheticSecondaries=0,tailDeltaCount=0;
  double tailDeltaEnergy=0;
};

class Gun final:public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{1};
  EventData* data;
public:
  explicit Gun(EventData* d):data(d){}
  void GeneratePrimaries(G4Event* event) override {
    *data=EventData{};
    data->event=event->GetEventID();
    G4Random::setTheSeed(static_cast<long>((cfg.seed+104729LL*data->event)%2147483646LL+1));
    const std::array<int,4> states{{0,1,15,27}};
    data->N=cfg.N<0?states[data->event%states.size()]:cfg.N;
    if(data->N>cfg.Z)throw std::runtime_error("Requested initial state N>Z");
    auto* ion=G4IonTable::GetIonTable()->GetIon(cfg.Z,cfg.A,0);
    data->initial=cfg.eu*ion->GetPDGMass()/kAtomicMassMeV;
    data->initialQ=cfg.Z-data->N;
    gun.SetParticleDefinition(ion);
    gun.SetParticleCharge(data->initialQ*eplus);
    gun.SetParticleEnergy(data->initial);
    gun.SetParticlePosition({0,0,-(.5*cfg.thickness+1)*um});
    gun.SetParticleMomentumDirection({0,0,1});
    gun.GeneratePrimaryVertex(event);
  }
};

class Events final:public G4UserEventAction {
  EventData* data;
public:
  explicit Events(EventData* d):data(d){}
  ~Events() override {delete data;}
  void EndOfEventAction(const G4Event* event) override {
    if(event->IsAborted())throw std::runtime_error("Aborted event");
    if(!std::isfinite(data->finalK)||!std::isfinite(data->edepAll)||data->edepAll<0||
       data->finalK<0||data->finalQ<0||data->finalQ>cfg.Z||data->maxDeltaRatio>1.000000001)
      throw std::runtime_error("Nonphysical event result");
    std::lock_guard<std::mutex> lock(resultsMutex);
    ++completed;
    maxPrimaryResidual=std::max(maxPrimaryResidual,std::abs(data->primaryResidual));
    output<<data->event<<','<<cfg.Z<<','<<cfg.A<<','<<data->N<<','<<cfg.eu<<','
      <<data->initial/MeV<<','<<data->initialQ<<','<<data->edepSlab/MeV<<','
      <<data->edepAll/MeV<<','<<(data->initial-data->finalK)/MeV<<','
      <<data->deltaCount<<','<<data->deltaEnergy/MeV<<','<<data->chargeExchange<<','
      <<data->scatter<<','<<data->finalK/MeV<<','<<data->finalQ<<','
      <<data->finalAngle<<','<<data->primarySteps<<','<<data->primaryResidual/MeV<<','
      <<data->maxStepResidual/MeV<<','<<data->escapedK/MeV<<','<<data->maxDeltaRatio<<','
      <<data->injectedK/MeV<<','<<data->syntheticSecondaries<<','<<data->secondaryInitialQ<<','
      <<data->secondaryFinalQ<<','<<data->secondaryFinalK/MeV<<','<<data->tailDeltaCount<<','
      <<data->tailDeltaEnergy/MeV<<'\n';
  }
};

class Steps final:public G4UserSteppingAction {
  EventData* data;
  G4ionEffectiveCharge expectedEffective;
public:
  explicit Steps(EventData* d):data(d){}
  void UserSteppingAction(const G4Step* step) override {
    const auto* track=step->GetTrack();
    const auto* pre=step->GetPreStepPoint();
    const auto* post=step->GetPostStepPoint();
    data->edepAll+=step->GetTotalEnergyDeposit();
    if(pre->GetPhysicalVolume()&&pre->GetPhysicalVolume()->GetName()=="Slab")
      data->edepSlab+=step->GetTotalEnergyDeposit();
    if(post->GetStepStatus()==fWorldBoundary)data->escapedK+=post->GetKineticEnergy();
    if(cfg.trace&&track->GetDefinition()->GetAtomicNumber()==cfg.Z&&
       track->GetDefinition()->GetAtomicMass()==cfg.A) {
      double chargeEnergy=post->GetKineticEnergy()+step->GetNonIonizingEnergyDeposit();
      for(const auto* child:*step->GetSecondaryInCurrentStep())
        if(child->GetDefinition()==G4Electron::Electron()&&child->GetCreatorProcess()&&
           child->GetCreatorProcess()->GetProcessName()=="CCIonDeltaElectron")
          chargeEnergy+=child->GetKineticEnergy();
      const double expected=expectedEffective.EffectiveCharge(track->GetDefinition(),
        pre->GetMaterial(),chargeEnergy)/eplus;
      const auto* defining=post->GetProcessDefinedStep();
      std::lock_guard<std::mutex> lock(resultsMutex);
      stepOutput<<data->event<<','<<track->GetTrackID()<<','<<track->GetParentID()<<','
        <<track->GetCurrentStepNumber()<<','<<G4CCIonChargeState::KineticEnergyPerNucleon(track->GetDefinition(),pre->GetKineticEnergy())/MeV<<','
        <<G4CCIonChargeState::KineticEnergyPerNucleon(track->GetDefinition(),post->GetKineticEnergy())/MeV<<','<<pre->GetCharge()/eplus<<','
        <<post->GetCharge()/eplus<<','<<track->GetDynamicParticle()->GetCharge()/eplus<<','
        <<expected<<','<<G4CCIonChargeState::KineticEnergyPerNucleon(track->GetDefinition(),chargeEnergy)/MeV<<','
        <<(defining?defining->GetProcessName():G4String())<<','
        <<pre->GetMaterial()->GetName()<<','<<step->GetStepLength()/mm<<'\n';
    }
    if(track->GetParentID()!=0)return;
    ++data->primarySteps;
    const auto* process=post->GetProcessDefinedStep();
    const auto name=process?process->GetProcessName():G4String();
    if(name=="CCIonChargeExchange")++data->chargeExchange;
    if(name=="CCTabulatedIonScattering"||name=="CoulombScat")++data->scatter;
    double secondaryK=0;
    const auto* secondaries=step->GetSecondaryInCurrentStep();
    const double M=track->GetParticleDefinition()->GetPDGMass();
    double bornDeltaK=0;
    for(const auto* secondary:*secondaries) {
      secondaryK+=secondary->GetKineticEnergy();
      if(secondary->GetParticleDefinition()==G4Electron::Electron()&&
         secondary->GetCreatorProcess()&&secondary->GetCreatorProcess()->GetProcessName()=="CCIonDeltaElectron") {
        ++data->deltaCount;data->deltaEnergy+=secondary->GetKineticEnergy();
        if(secondary->GetKineticEnergy()<1*keV) {
          ++data->tailDeltaCount;data->tailDeltaEnergy+=secondary->GetKineticEnergy();
        }
        bornDeltaK+=secondary->GetKineticEnergy();
      }
    }
    if(bornDeltaK>0) {
      const double K=post->GetKineticEnergy()+bornDeltaK;
      const double tau=K/M,r=electron_mass_c2/M;
      const double tmax=2*electron_mass_c2*tau*(tau+2)/(1+2*(tau+1)*r+r*r);
      data->maxDeltaRatio=std::max(data->maxDeltaRatio,bornDeltaK/tmax);
    }
    const double residual=pre->GetKineticEnergy()-post->GetKineticEnergy()-
      step->GetTotalEnergyDeposit()-secondaryK;
    if(std::abs(residual)>1e-5*MeV) {
      std::lock_guard<std::mutex> lock(resultsMutex);
      std::cerr<<std::setprecision(17)<<"PRIMARY_STEP_IMBALANCE event="<<data->event
        <<" step="<<track->GetCurrentStepNumber()<<" process="<<name
        <<" pre_MeV="<<pre->GetKineticEnergy()/MeV<<" post_MeV="<<post->GetKineticEnergy()/MeV
        <<" edep_MeV="<<step->GetTotalEnergyDeposit()/MeV
        <<" niel_MeV="<<step->GetNonIonizingEnergyDeposit()/MeV<<" secondary_MeV="<<secondaryK/MeV
        <<" residual_MeV="<<residual/MeV<<" length_mm="<<step->GetStepLength()/mm
        <<" z_mm="<<post->GetPosition().z()/mm<<" pre_q="<<pre->GetCharge()/eplus
        <<" post_q="<<post->GetCharge()/eplus<<" couple="<<pre->GetMaterialCutsCouple()->GetIndex()
        <<" material="<<pre->GetMaterial()->GetName()<<'\n';
    }
    data->primaryResidual+=residual;
    data->maxStepResidual=std::max(data->maxStepResidual,std::abs(residual));
  }
};

class Tracking final:public G4UserTrackingAction {
  EventData* data;
public:
  explicit Tracking(EventData* d):data(d){}
  void PostUserTrackingAction(const G4Track* track) override {
    if(track->GetParentID()!=0) {
      if(cfg.secondary&&track->GetDefinition()->GetAtomicNumber()==cfg.Z&&
         track->GetDefinition()->GetAtomicMass()==cfg.A) {
        data->secondaryFinalK=track->GetKineticEnergy();
        data->secondaryFinalQ=track->GetDynamicParticle()->GetCharge()/eplus;
      }
      return;
    }
    data->finalK=track->GetKineticEnergy();
    data->finalQ=track->GetDynamicParticle()->GetCharge()/eplus;
    data->finalAngle=std::acos(std::clamp(track->GetMomentumDirection().z(),-1.,1.));
    if(cfg.secondary) {
      const std::array<int,4> states{{0,1,15,27}};
      const int n=states[(data->event+1)%states.size()];
      auto* ion=G4IonTable::GetIonTable()->GetIon(cfg.Z,cfg.A,0);
      const double energy=(cfg.secondaryEu>0?cfg.secondaryEu:cfg.eu)*ion->GetPDGMass()/kAtomicMassMeV;
      auto* dynamic=new G4DynamicParticle(ion,{0,0,1},energy);
      dynamic->SetCharge((cfg.Z-n)*eplus);
      auto* child=new G4Track(dynamic,track->GetGlobalTime(),{0,0,-(.5*cfg.thickness+1)*um});
      child->SetParentID(track->GetTrackID());
      fpTrackingManager->GimmeSecondaries()->push_back(child);
      ++data->syntheticSecondaries;
      data->injectedK+=energy;
      data->secondaryInitialQ=cfg.Z-n;
    }
  }
};

class Actions final:public G4VUserActionInitialization {
  void Build() const override {
    auto* data=new EventData;
    SetUserAction(new Gun(data));
    SetUserAction(new Events(data));
    SetUserAction(new Steps(data));
    SetUserAction(new Tracking(data));
  }
};

void Audit() {
  std::ofstream audit(cfg.out+".audit.txt");
  audit<<std::setprecision(17);
  auto* pm=G4GenericIon::GenericIon()->GetProcessManager();
  std::map<std::string,int> names;
  for(int i=0;i<pm->GetProcessListLength();++i) {
    auto* p=(*pm->GetProcessList())[i];
    if(pm->GetProcessActivation(p)) {
      ++names[p->GetProcessName()];
      audit<<"ALONG_ORDER "<<p->GetProcessName()<<' '<<pm->GetProcessOrdering(p,idxAlongStep)<<'\n';
    }
  }
  for(const auto& entry:names)audit<<"PROCESS "<<entry.first<<' '<<entry.second<<'\n';
  for(const auto* name:{"CCIonIoni","CCIonDeltaElectron","CCIonChargeExchange","CCNonEMChargeReset"})
    if(names[name]!=1)throw std::runtime_error(std::string("Expected exactly one ")+name);
  for(const auto* name:{"ionIoni","ionInelastic","hadElastic"})
    if(names[name])throw std::runtime_error(std::string("Unexpected process ")+name);
  auto* cuts=G4ProductionCutsTable::GetProductionCutsTable();
  for(std::size_t i=0;i<cuts->GetTableSize();++i) {
    auto* couple=cuts->GetMaterialCutsCouple(i);
    if(couple->GetMaterial()!=slabMaterial)continue;
    const auto energyCut=cuts->GetEnergyCutsVector(1)->at(i)/keV;
    audit<<"ELECTRON_CUT_KEV "<<energyCut<<'\n';
    if(cfg.cutKeV>0&&std::abs(energyCut-cfg.cutKeV)>1e-12)
      throw std::runtime_error("Requested exact electron cut was not applied");
  }
  for(double e:{29.,30.,35.,40.,60.,80.,100.,400.,700.,1000.,1200.,1500.,1750.,2000.,2100.})
    audit<<"WEIGHT "<<e<<' '<<G4CCCurrentChargeConfig::GlobalTransitionWeight(e)<<' '
      <<G4CCCurrentChargeConfig::IonFluctuationMixProbability(e)<<'\n';
}

void KernelAudit() {
  if(!cfg.kernel||cfg.legacy)return;
  auto* ion=G4IonTable::GetIonTable()->GetIon(cfg.Z,cfg.A,0);
  if(!G4CCDiracPhysics::Applies(ion,slabMaterial))return;
  G4CCLindhardSorensenIonModel* model=nullptr;
  auto* pm=G4GenericIon::GenericIon()->GetProcessManager();
  for(int i=0;i<pm->GetProcessListLength();++i) {
    auto* loss=dynamic_cast<G4VEnergyLossProcess*>((*pm->GetProcessList())[i]);
    if(!loss||loss->GetProcessName()!="CCIonIoni")continue;
    for(std::size_t j=0;j<loss->NumberOfModels();++j)
      if(auto* candidate=dynamic_cast<G4CCLindhardSorensenIonModel*>(loss->GetModelByIndex(j)))model=candidate;
  }
  if(!model)throw std::runtime_error("Cannot find initialized CC LS model for kernel audit");
  auto* cuts=G4ProductionCutsTable::GetProductionCutsTable();
  for(std::size_t i=0;i<cuts->GetTableSize();++i) {
    auto* couple=cuts->GetMaterialCutsCouple(i);
    if(couple->GetMaterial()==slabMaterial){model->SetCurrentCouple(couple);break;}
  }
  std::vector<double> energies{20.,80.,100.,400.,1000.,2500.};
  for(const double e:cfg.edges)for(const double scale:{1.-1e-8,1.,1.+1e-8})energies.push_back(e*scale);
  std::sort(energies.begin(),energies.end());
  energies.erase(std::unique(energies.begin(),energies.end()),energies.end());
  std::ofstream probe(cfg.out+".kernel.csv");
  probe<<std::setprecision(17)<<"E_u,N,cut_keV,weight,ion_fluct_probability,full_MeV_mm,hard_MeV_mm,restricted_MeV_mm,balance_relative,dirac_full_MeV_mm,dirac_hard_MeV_mm,GLOBAL_total_cm2,fallback_form_factor_half_Tmax\n";
  for(double e:energies)for(int n:{0,1,15,27})for(double cut:{.99*keV,1*keV,10*keV}) {
    if(n>cfg.Z)continue;
    const double K=e*ion->GetPDGMass()/kAtomicMassMeV;
    double restricted=0,full=0,hard=0;
    if(!model->ComputeDiracRestrictedDEDX(slabMaterial,ion,K,cut,cfg.Z-n,restricted,&full,&hard))
      throw std::runtime_error("Dirac restricted model declined covered kernel case");
    const double w=G4CCCurrentChargeConfig::GlobalTransitionWeight(e);
    const double error=std::abs(restricted+hard-full)/std::max(full,1e-300);
    if(error>2e-13||restricted<0||!std::isfinite(restricted))
      throw std::runtime_error("Kernel stopping decomposition failed");
    double dfull=0,dhard=0;
    if(e>=30&&e<=2000) {
      dfull=G4CCDiracPhysics::FullStopping(slabMaterial,ion,n,K);
      dhard=G4CCDiracPhysics::HardStopping(slabMaterial,ion,n,K,cut);
      if(w==1 && (std::abs(dfull-full)>1e-12*full||std::abs(dhard-hard)>1e-12*std::max(hard,1e-300)))
        throw std::runtime_error("Plateau changed the direct tabulated spectrum or full stopping");
    }
    double global=0;
    const auto tr=G4CCGlobalChargeExchangeTable::Instance().FindTransitions(cfg.Z,cfg.A,14,e,n);
    for(const auto& t:tr)global+=t.sigmaCm2;
    const double tau=K/ion->GetPDGMass(),ratio=electron_mass_c2/ion->GetPDGMass();
    const double top=2*electron_mass_c2*tau*(tau+2)/(1+2*(tau+1)*ratio+ratio*ratio);
    const double ff=G4CCIonDeltaElectronProcess::FormFactorAcceptance(ion,K,.5*top);
    probe<<e<<','<<n<<','<<cut/keV<<','<<w<<','<<G4CCCurrentChargeConfig::IonFluctuationMixProbability(e)<<','
      <<full/(MeV/mm)<<','<<hard/(MeV/mm)<<','<<restricted/(MeV/mm)<<','<<error<<','
      <<dfull/(MeV/mm)<<','<<dhard/(MeV/mm)<<','<<global<<','<<ff<<'\n';
  }
  std::ofstream samples(cfg.out+".physical_samples.csv");
  samples<<std::setprecision(17)<<"E_u,N,cut_keV,samples,expected_mean_MeV,sample_mean_MeV,mean_standard_errors,max_T_over_Tmax\n";
  std::mt19937_64 rng(20261101);std::uniform_real_distribution<double> random(0,1);
  for(double e:{35.,400.,1750.})for(int n:{0,1,15,27})for(double cut:{.99*keV,1*keV,10*keV}) {
    if(n>cfg.Z)continue;
    const double K=e*ion->GetPDGMass()/kAtomicMassMeV;
    const double tau=K/ion->GetPDGMass(),r=electron_mass_c2/ion->GetPDGMass();
    const double top=2*electron_mass_c2*tau*(tau+2)/(1+2*(tau+1)*r+r*r);
    const auto moments=G4CCDiracPhysics::Moments(ion,n,K,cut);
    const double expected=moments[1]/moments[0];
    const double variance=std::max(moments[2]/moments[0]-expected*expected,0.);
    constexpr int count=16384;
    double sum=0,maximum=0;
    for(int i=0;i<count;++i) {
      const double t=G4CCDiracPhysics::Sample(ion,n,K,cut,random(rng),random(rng));
      if(t<cut*(1-1e-12)||t>top*(1+1e-12))throw std::runtime_error("Finite-mass sample out of bounds");
      sum+=t/MeV;maximum=std::max(maximum,t/top);
    }
    const double score=std::abs(sum/count-expected)/std::sqrt(variance/count);
    if(score>6)throw std::runtime_error("Finite-mass sampled mean failed");
    samples<<e<<','<<n<<','<<cut/keV<<','<<count<<','<<expected<<','<<sum/count<<','<<score<<','<<maximum<<'\n';
  }
}
}

int main(int argc,char** argv) {
  try {
    Parse(argc,argv);
    output.open(cfg.out);
    if(!output)throw std::runtime_error("Cannot open output CSV");
    output<<std::setprecision(17)<<"event_id,Z,A,N_initial,E_u_MeV,K_initial_MeV,q_initial,edep_slab_MeV,edep_all_MeV,primary_loss_MeV,delta_count,delta_energy_MeV,charge_exchange_count,ss_count,K_final_MeV,q_final,angle_final_rad,primary_steps,primary_balance_MeV,max_step_balance_MeV,escaped_kinetic_MeV,max_delta_T_over_Tmax,injected_kinetic_MeV,synthetic_secondary_count,secondary_q_initial,secondary_q_final,secondary_K_final_MeV,tail_delta_count,tail_delta_energy_MeV\n";
    if(cfg.trace) {
      stepOutput.open(cfg.out+".steps.csv");
      if(!stepOutput)throw std::runtime_error("Cannot open step trace");
      stepOutput<<std::setprecision(17)<<"event_id,track_id,parent_id,step,E_pre_u,E_post_u,q_pre,q_post,q_track,q_effective,E_charge_u,process,material,length_mm\n";
    }
    TestExceptionHandler handler;
    std::unique_ptr<G4RunManager> manager;
    if(cfg.threads>1) {
      auto* mt=new G4MTRunManager;
      mt->SetNumberOfThreads(cfg.threads);
      manager.reset(mt);
    } else manager=std::make_unique<G4RunManager>();
    G4StateManager::GetStateManager()->SetExceptionHandler(&handler);
    manager->SetVerboseLevel(0);manager->SetPrintProgress(0);
    manager->SetUserInitialization(new Detector);
    manager->SetUserInitialization(new Physics);
    manager->SetUserInitialization(new Actions);
    manager->Initialize();
    manager->BeamOn(cfg.events);
    Audit();
    KernelAudit();
    if(completed!=static_cast<std::size_t>(cfg.events))throw std::runtime_error("Wrong completed event count");
    output.close();
    double peakMB=0;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS memory{};
    if(GetProcessMemoryInfo(GetCurrentProcess(),&memory,sizeof(memory)))peakMB=memory.PeakWorkingSetSize/(1024.*1024.);
#endif
    std::cout<<"DIRAC_TRANSPORT_PASS events="<<completed<<" threads="<<cfg.threads
      <<" wall_s="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
      <<" peak_working_set_MB="<<peakMB<<" max_primary_balance_MeV="<<maxPrimaryResidual/MeV<<'\n';
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"DIRAC_TRANSPORT_ERROR "<<e.what()<<'\n';return 1;
  }
}
