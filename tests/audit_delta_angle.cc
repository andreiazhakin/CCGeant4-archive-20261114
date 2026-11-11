#define main retained_transport_main
#include "dirac_transport.cc"
#undef main

int main(int argc,char** argv) {
  try {
    Parse(argc,argv);
    TestExceptionHandler exceptions;
    G4StateManager::GetStateManager()->SetExceptionHandler(&exceptions);
    auto manager=std::make_unique<G4RunManager>();
    manager->SetVerboseLevel(0);manager->SetPrintProgress(0);
    manager->SetUserInitialization(new Detector);
    manager->SetUserInitialization(new Physics);
    manager->SetUserInitialization(new Actions);
    manager->Initialize();manager->BeamOn(0);
    auto* ion=G4IonTable::GetIonTable()->GetIon(54,124,0);
    G4CCIonDeltaElectronProcess* process=nullptr;
    const auto* processes=G4GenericIon::GenericIon()->GetProcessManager()->GetProcessList();
    for(std::size_t i=0;i<processes->size();++i)
      if(auto* p=dynamic_cast<G4CCIonDeltaElectronProcess*>((*processes)[i]))process=p;
    if(!process)throw std::runtime_error("No registered CC delta process");
    const G4MaterialCutsCouple* couple=nullptr;
    auto* cuts=G4ProductionCutsTable::GetProductionCutsTable();
    for(std::size_t i=0;i<cuts->GetTableSize();++i)
      if(cuts->GetMaterialCutsCouple(i)->GetMaterial()==slabMaterial)couple=cuts->GetMaterialCutsCouple(i);
    if(!couple)throw std::runtime_error("No silicon couple");
    std::ofstream file(cfg.out);
    if(!file)throw std::runtime_error("Cannot open output");
    file<<std::setprecision(17)<<"Eu,N,sample,count,T_MeV,dx,dy,dz,next_random,mean_free_path_mm\n";
    std::size_t generated=0;
    for(double eu:{35.,400.,1750.,2500.})for(int n:{1,15}) {
      G4Random::setTheSeed(1234567);
      G4Step step;
      auto* dynamic=new G4DynamicParticle(ion,{0,0,1},eu*MeV*ion->GetPDGMass()/(931.494*MeV));
      dynamic->SetCharge((54-n)*eplus);
      G4Track track(dynamic,0,{});track.SetTrackID(1);track.SetStep(&step);
      step.SetTrack(&track);
      step.GetPreStepPoint()->SetMaterial(slabMaterial);
      step.GetPreStepPoint()->SetMaterialCutsCouple(couple);
      process->StartTracking(&track);
      G4ForceCondition condition=NotForced;
      const double meanFreePath=process->GetMeanFreePath(track,0.,&condition);
      if(!std::isfinite(meanFreePath)||meanFreePath<=0.)throw std::runtime_error("Invalid mean free path");
      for(int sample=0;sample<128;++sample) {
        auto* change=process->PostStepDoIt(track,step);
        const int count=change->GetNumberOfSecondaries();
        if(count<0||count>1)throw std::runtime_error("Unexpected secondary count");
        double energy=0;G4ThreeVector direction;
        if(count) {
          const auto* electron=change->GetSecondary(0);
          if(electron->GetDefinition()!=G4Electron::Electron())throw std::runtime_error("Non-electron secondary");
          energy=electron->GetKineticEnergy();direction=electron->GetMomentumDirection();
          if(!std::isfinite(energy)||energy<=0||std::abs(direction.mag2()-1)>1e-12)
            throw std::runtime_error("Invalid secondary kinematics");
          ++generated;delete electron;
        }
        change->Clear();
        file<<eu<<','<<n<<','<<sample<<','<<count<<','<<energy/MeV<<','
            <<direction.x()<<','<<direction.y()<<','<<direction.z()<<','<<G4UniformRand()<<','<<meanFreePath/mm<<'\n';
      }
      process->EndTracking();
    }
    if(!generated)throw std::runtime_error("No deltas sampled");
    std::cout<<"AUDIT_DELTA_ANGLE_PASS samples=1024 generated="<<generated<<'\n';
    return 0;
  }catch(const std::exception& e){std::cerr<<"AUDIT_DELTA_ANGLE_FAIL "<<e.what()<<'\n';return 1;}
}
