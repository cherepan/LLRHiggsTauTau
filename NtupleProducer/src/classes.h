//#include <AnalysisDataFormats/CMGTools/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include <vector>
#include <DataFormats/PatCandidates/interface/PFParticle.h>

edm::Ptr<pat::PFParticle> dummy1;
pat::UserHolder<std::vector<edm::Ptr<pat::PFParticle> > > dummy2;
pat::UserHolder<std::vector<double, std::allocator<double> > > dummy3;
pat::UserHolder<std::vector<std::vector<double > > > dummy4;
std::vector< std::vector<std::vector<double> > >    dummy5;
std::vector< std::vector<int > >    dummy6;
std::vector< std::vector<unsigned int> >     dummy7;
