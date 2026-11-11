#define main retained_transport_main
#include "dirac_transport.cc"
#undef main
#include "G4CCAtimaEquilibriumDedxTable.hh"
#include <tuple>

int main(int argc,char** argv) {
  try {
    Parse(argc,argv);
    TestExceptionHandler handler;
    G4StateManager::GetStateManager()->SetExceptionHandler(&handler);
    auto manager=std::make_unique<G4RunManager>();
    manager->SetVerboseLevel(0);manager->SetPrintProgress(0);
    manager->SetUserInitialization(new Detector);
    manager->SetUserInitialization(new Physics);
    manager->SetUserInitialization(new Actions);
    manager->Initialize();
    auto* ion=G4IonTable::GetIonTable()->GetIon(54,124,0);
    const auto& low=G4CCAtimaEquilibriumDedxTable::Instance();
    struct Entry {G4Material* material;int z;};
    std::vector<Entry> entries;
    for(const auto& item:std::vector<std::tuple<std::string,int,double>>{
        {"ReferenceH",1,low.GetReferenceAtomicMass(54,1)},
        {"Deuterium",1,2.01410177812},{"Tritium",1,3.01604928199},
        {"ReferenceSi",14,low.GetReferenceAtomicMass(54,14)},
        {"Silicon28",14,27.97692653465},{"Silicon29",14,28.97649466490}}) {
      entries.push_back({new G4Material(std::get<0>(item),std::get<1>(item),
        std::get<2>(item)*g/mole,1*g/cm3),std::get<1>(item)});
    }
    entries.push_back({G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"),14});
    std::ofstream out(cfg.out);
    out<<std::setprecision(17)<<"material,Zt,Eu,atomic_density_mm3,low_MeV_mm,low_per_atom_MeV_mm2,upper_MeV_mm,upper_over_low\n";
    std::map<int,double> boundaryRatio;
    double maxPerAtomError=0,maxBoundarySpread=0,maxBoundaryOffset=0;
    for(const auto& item:entries)for(double eu:{2.,10.,29.999,30.,40.}) {
      auto* material=item.material;
      double k=eu/931.494*ion->GetPDGMass();
      while((eu==2.||eu==30.) && G4CCIonChargeState::KineticEnergyPerNucleon(ion,k)/MeV<eu)
        k=std::nextafter(k,std::numeric_limits<double>::infinity());
      while(eu==40. && G4CCIonChargeState::KineticEnergyPerNucleon(ion,k)/MeV>eu)
        k=std::nextafter(k,0.);
      double value=0,mass=0;
      if(!low.GetMaterialDEDXPerVolume(ion,material,k,value)||
         !low.GetMassDEDX(54,item.z,G4CCIonChargeState::KineticEnergyPerNucleon(ion,k)/MeV,mass))
        throw std::runtime_error("Low-energy pair declined");
      const double density=material->GetVecNbOfAtomsPerVolume()[0];
      const double perAtom=value/density;
      const double expected=mass*low.GetReferenceAtomicMass(54,item.z)/(Avogadro*mole)*MeV*cm2;
      maxPerAtomError=std::max(maxPerAtomError,std::abs(perAtom/expected-1));
      double upper=0,ratio=0;
      if(eu>=30.) {
        upper=G4CCDiracPhysics::EquilibriumStopping(material,ion,k);ratio=upper/value;
        if(eu==30.) {
          if(!boundaryRatio.count(item.z))boundaryRatio[item.z]=ratio;
          maxBoundarySpread=std::max(maxBoundarySpread,std::abs(ratio-boundaryRatio[item.z]));
          maxBoundaryOffset=std::max(maxBoundaryOffset,std::abs(ratio-1));
        }
      }
      out<<material->GetName()<<','<<item.z<<','<<eu<<','<<density*mm3<<','<<value/(MeV/mm)<<','
         <<perAtom/(MeV*mm2)<<','<<upper/(MeV/mm)<<','<<ratio<<'\n';
    }
    auto* water=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    double full=0,sum=0;const double k=10./931.494*ion->GetPDGMass();
    low.GetMaterialDEDXPerVolume(ion,water,k,full);
    for(std::size_t i=0;i<water->GetNumberOfElements();++i) {
      const int z=(*water->GetElementVector())[i]->GetZasInt();double v=0;
      low.GetMassDEDX(54,z,10.,v);
      sum+=water->GetVecNbOfAtomsPerVolume()[i]*v*low.GetReferenceAtomicMass(54,z)/(Avogadro*mole)*MeV*cm2;
    }
    const double mixtureError=std::abs(full/sum-1);
    if(maxPerAtomError>1e-12||maxBoundarySpread>1e-12||maxBoundaryOffset>1e-4||mixtureError>1e-12)
      throw std::runtime_error("Isotopic per-atom normalization/boundary failure");
    std::cout<<std::setprecision(17)<<"ISOTOPE_PASS max_per_atom_error="<<maxPerAtomError
      <<" max_boundary_isotope_spread="<<maxBoundarySpread<<" max_boundary_offset="<<maxBoundaryOffset
      <<" mixture_error="<<mixtureError<<std::endl;
    return 0;
  }catch(const std::exception& e){std::cerr<<"ISOTOPE_ERROR "<<e.what()<<std::endl;return 1;}
}
