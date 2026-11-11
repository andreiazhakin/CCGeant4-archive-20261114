#include "G4CCDiracTable.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

int main(int argc,char** argv) {
  try {
    if(argc!=5)throw std::runtime_error("Usage: table_kernel_probe DELTA STOPPING Z OUTPUT.csv");
    auto& table=G4CCDiracTable::Instance();
    table.Configure(argv[1],argv[2]);table.Validate();
    const int z=std::stoi(argv[3]);
    std::ofstream out(argv[4]);
    if(!out)throw std::runtime_error("Cannot write probe CSV");
    out<<std::setprecision(17)<<"Z,Zt,N,E,cut,M0,M1,M2,full,eq,B,samples,mean_sample,mean_sample2,min_sample,max_sample\n";
    std::mt19937_64 rng(11012026);
    std::uniform_real_distribution<double> random(0,1);
    for(double energy:{30.,30.25,35.,40.,100.,123.45,400.,699.5,1500.,1999.5,2000.}) {
      const double tmax=table.MaxTransfer(energy);
      for(int n:{0,1,15,27})for(double cut:{.001,.01,.1,.9*tmax}) {
        const auto moments=table.DeltaMoments(z,n,energy,cut);
        double sum=0,sum2=0,minimum=tmax,maximum=0;
        int samples=(energy==30.25||energy==400.||energy==1999.5)?16384:0;
        if(cut>=tmax)samples=0;
        for(int i=0;i<samples;++i) {
          const double t=table.SampleDelta(z,n,energy,cut,random(rng),random(rng));
          if(t<cut*(1-1e-12)||t>tmax*(1+1e-12))throw std::runtime_error("Sample outside spectrum support");
          sum+=t;sum2+=t*t;minimum=std::min(minimum,t);maximum=std::max(maximum,t);
        }
        out<<z<<",14,"<<n<<','<<energy<<','<<cut<<','<<moments[0]<<','<<moments[1]<<','<<moments[2]<<','
          <<table.FullStopping(z,14,n,energy)<<','<<table.EquilibriumStopping(z,14,energy)<<','
          <<table.SoftCoefficient(z,14,energy)<<','<<samples<<','<<(samples?sum/samples:0)<<','
          <<(samples?sum2/samples:0)<<','<<(samples?minimum:0)<<','<<maximum<<'\n';
      }
    }
    std::cout<<"TABLE_KERNEL_PASS memory_bytes="<<table.MemoryBytes()<<'\n';
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
