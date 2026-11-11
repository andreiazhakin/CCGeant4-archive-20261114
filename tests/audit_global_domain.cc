#include "G4CCCurrentChargeConfig.hh"
#include "G4CCPhysicsOptions.hh"
#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCIonChargeState.hh"
#include "G4CCIonChargeExchangeProcess.hh"
#include "G4CCNonEMChargeResetProcess.hh"
#include "G4CCDiracPhysics.hh"
#include "G4DynamicParticle.hh"
#include "G4Ions.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VExceptionHandler.hh"
#include "G4VDiscreteProcess.hh"
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

class Exceptions final:public G4VExceptionHandler {
  G4bool Notify(const char* where,const char* code,G4ExceptionSeverity severity,const char* message) override {
    std::cerr<<"G4_EXCEPTION "<<where<<' '<<code<<' '<<message<<'\n';
    if(severity==FatalException||severity==FatalErrorInArgument)throw std::runtime_error(code);
    return false;
  }
};
class FakeReaction final:public G4VDiscreteProcess {
 public:
  FakeReaction():G4VDiscreteProcess("AuditReaction",fHadronic){}
  G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*) override{return DBL_MAX;}
};
int main(int argc,char** argv) {
  try {
    if(argc!=8)throw std::runtime_error("MODE GLOBAL DELTA STOPPING LOW Zp Zt required");
    Exceptions errors;G4StateManager::GetStateManager()->SetExceptionHandler(&errors);
    G4CCPhysicsOptions options;
    options.tables.global=argv[2];options.tables.delta=argv[3];
    options.tables.stopping=argv[4];options.tables.atimaLow=argv[5];
    G4CCCurrentChargeConfig::Configure(options);
    const int z=std::stoi(argv[6]),zt=std::stoi(argv[7]);
    auto& global=G4CCGlobalChargeExchangeTable::Instance();
    if(std::string(argv[1])=="table") {
      global.EnsureLoaded();global.RequirePair(z,zt);
      double sum=0;
      for(const auto& tr:global.FindTransitions(z,124,zt,400.,0))sum+=tr.sigmaCm2;
      std::cout<<std::setprecision(17)<<"TABLE_PASS curves="<<global.NumberOfCurves()
        <<" points="<<global.NumberOfPoints()<<" sigma_N0_E400="<<sum<<'\n';
      for(const double energy:{30.0001,40.,400.,1500.,1999.999})for(int n=0;n<=27;++n)
        for(const auto& tr:global.FindTransitions(z,124,zt,energy,n))
          std::cout<<"RATE "<<energy<<' '<<n<<' '<<tr.dQe<<' '<<tr.sigmaCm2<<'\n';
      return 0;
    }
    const int a=z==20?40:(z==29?63:(z==54?124:237));
    auto* ion=new G4Ions("audit_ion",a*931.494*MeV,0.,z*eplus,0,1,0,0,0,0,
      "nucleus",0,a,1000000000+z*10000+a*10,true,-1.,nullptr,false,"generic");
    auto* material=zt==14?G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"):
      new G4Material("AuditTarget",static_cast<double>(zt),std::max(1.,2.5*zt)*g/mole,1*g/cm3);
    const bool ce=G4CCIonChargeState::CurrentChargeModelApplies(ion,material);
    const bool em=G4CCIonChargeState::IonisationModelApplies(ion,material);
    if(ce!=em)throw std::runtime_error("Explicit domains disagree");
    auto* dynamic=new G4DynamicParticle(ion,{0,0,1},400.*ion->GetPDGMass()/931.494);
    dynamic->SetCharge((z-1)*eplus);
    G4Track track(dynamic,0.,{});G4Step step;
    step.SetTrack(&track);track.SetStep(&step);
    step.GetPreStepPoint()->SetMaterial(material);
    step.GetPostStepPoint()->SetMaterial(material);
    FakeReaction reaction;step.GetPostStepPoint()->SetProcessDefinedStep(&reaction);
    G4CCIonChargeExchangeProcess exchange;
    G4ForceCondition condition=NotForced;
    const double mfp=exchange.GetMeanFreePath(track,0.,&condition);
    G4CCNonEMChargeResetProcess reset;
    auto* change=static_cast<G4ParticleChangeForLoss*>(reset.PostStepDoIt(track,step));
    const double q=change->GetProposedCharge()/eplus;
    if(std::abs(q-(ce?z:z-1))>1e-12)throw std::runtime_error("Reset domain disagrees");
    if((ce&&!(mfp>0.&&mfp<DBL_MAX))||(!ce&&mfp!=DBL_MAX))throw std::runtime_error("CE domain disagrees");
    std::cout<<std::setprecision(17)<<"DOMAIN_PASS CE="<<ce<<" EM="<<em<<" reset_q="<<q
      <<" mfp_mm="<<mfp/mm<<" max_target="<<G4CCIonChargeState::MaxTargetZ()<<'\n';
    return 0;
  }catch(const std::exception& e){std::cerr<<"AUDIT_ERROR "<<e.what()<<'\n';return 1;}
}
