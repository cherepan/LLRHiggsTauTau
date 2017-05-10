/** \class TauFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author G. Ortona (LLR)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "BDTId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include <vector>
#include <string>
#include <TMath.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;

class TauFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller(){
    
  };  
  //ByIsolationMVA3oldDMwoLTraw
 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<pat::TauRefVector> theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  edm::EDGetTokenT<vector<Vertex> > theVtxTag ;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> thePFCandTag;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotTag;

  const std::string theDiscriminatorTag;
  const StringCutObjectSelector<pat::Tau, true> cut;
  const CutSet<pat::Tau> flags;
  // const std::string NominalUpOrDown;
  const double NominalTESCorrection; // value of correction of centrale TES value
  const bool ApplyTESCentralCorr; // shift the central TES value
  // const bool ApplyTESUpDown; // compute Up/Down TES variation

  vector<string> tauIntDiscrims_; // tau discrims to be added as userInt
  vector<string> tauFloatDiscrims_; // tau discrims to be added as userFloats
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(consumes<pat::TauRefVector>(iConfig.getParameter<InputTag>("src"))),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  thePFCandTag(consumes<edm::View<pat::PackedCandidate>>        (iConfig.getParameter<edm::InputTag>("PFCollection"))),
  beamSpotTag(consumes<reco::BeamSpot>                         (iConfig.getParameter<edm::InputTag>("offlinebeamSpot"))),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")), 

  // NominalUpOrDown(iConfig.getParameter<std::string>("NominalUpOrDown")),
  NominalTESCorrection(iConfig.getParameter<double>("NominalTESCorrection")),
  ApplyTESCentralCorr(iConfig.getParameter<bool>("ApplyTESCentralCorr"))
  // ApplyTESUpDown(iConfig.getParameter<bool>("ApplyTESUpDown"))
{
  produces<pat::TauCollection>();

  tauIntDiscrims_ = 
  {
    "decayModeFinding", // it is decayModeFindingOldDMs
    "decayModeFindingNewDMs",
    
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    
    "byVLooseIsolationMVArun2v1DBoldDMwLT",
    "byLooseIsolationMVArun2v1DBoldDMwLT",
    "byMediumIsolationMVArun2v1DBoldDMwLT",
    "byTightIsolationMVArun2v1DBoldDMwLT",
    "byVTightIsolationMVArun2v1DBoldDMwLT",

    "byVLooseIsolationMVArun2v1DBnewDMwLT",    
    "byLooseIsolationMVArun2v1DBnewDMwLT",
    "byMediumIsolationMVArun2v1DBnewDMwLT",
    "byTightIsolationMVArun2v1DBnewDMwLT",
    "byVTightIsolationMVArun2v1DBnewDMwLT",

    "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
    "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
    "byTightIsolationMVArun2v1DBdR03oldDMwLT",
    "byVTightIsolationMVArun2v1DBdR03oldDMwLT",
    
    "byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
    "byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",
    "byTightCombinedIsolationDeltaBetaCorr3HitsdR03",

    "againstElectronMVA5category",
    
    "byLooseIsolationMVA3newDMwLT",
    "byLooseIsolationMVA3oldDMwLT",
    "byLoosePileupWeightedIsolation3Hits",
    "byMediumIsolationMVA3newDMwLT",
    "byMediumIsolationMVA3oldDMwLT",
    "byMediumPileupWeightedIsolation3Hits",
    "byTightIsolationMVA3newDMwLT",
    "byTightIsolationMVA3oldDMwLT",
    "byTightPileupWeightedIsolation3Hits",
    
    "byVLooseIsolationMVA3newDMwLT",
    "byVTightIsolationMVA3newDMwLT",
    "byVVTightIsolationMVA3newDMwLT",

    "byVLooseIsolationMVA3oldDMwLT",
    "byVTightIsolationMVA3oldDMwLT",
    "byVVTightIsolationMVA3oldDMwLT",

    "againstMuonLoose3",
    "againstMuonTight3",

    "againstElectronVLooseMVA6",
    "againstElectronLooseMVA6",
    "againstElectronMediumMVA6",
    "againstElectronTightMVA6",
    "againstElectronVTightMVA6"

  };

  tauFloatDiscrims_ =
  {
    "byCombinedIsolationDeltaBetaCorrRaw3Hits",
    "byIsolationMVArun2v1DBoldDMwLTraw",
    "byIsolationMVA3oldDMwoLTraw",
    "byIsolationMVA3oldDMwLTraw",
    "byIsolationMVA3newDMwoLTraw",
    "againstElectronMVA5raw",
    "byPhotonPtSumOutsideSignalCone",
    "byPileupWeightedIsolationRaw3Hits",
    "footprintCorrection",
    "neutralIsoPtSumWeight",
    "photonPtSumOutsideSignalCone",
    "byIsolationMVA3newDMwLTraw",
    "chargedIsoPtSum",
    "neutralIsoPtSum",
    "puCorrPtSum",
  };


}


void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
//read one PFTauDiscriminator (set discriminatorSrc_ in to an edm::InputTag before)

  // Get leptons and discriminators
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByToken(theCandidateTag, tauHandle);
    
  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByToken(theVtxTag, vertexs);


  edm::Handle<edm::View<pat::PackedCandidate> >pfCandHandle;
  iEvent.getByToken(thePFCandTag,pfCandHandle);
  const edm::View<pat::PackedCandidate>* cands = pfCandHandle.product();

  // Output collection
  auto_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int itau = 0; itau < tauHandle->size(); ++itau){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[itau].get()));
    std::vector<double > SVPos;     
    std::vector<double > SVCov;     
    std::vector<std::vector<double> > iPionP4;
    std::vector<double> iPionCharge;
    std::vector<double> SVChi2NDofMatchingQual;

    // Nominal TES Correction
    double Shift = 1.+NominalTESCorrection/100.;
    double shiftP = 1.;
    double shiftMass = 1.;
    
    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. && ApplyTESCentralCorr)
    {
      if(l.decayMode()>=1 && l.decayMode()<=10){
        shiftP = Shift;
        shiftMass = Shift;
      }
      else if(l.decayMode()==0){
        shiftP = Shift;
        shiftMass = 1.;
      }
      // else if(l.decayMode()==10){
      //   shiftP = Shift;
      //   shiftMass = Shift;
      // }
    }
    
    double pxS_Nominal = l.px()*shiftP;
    double pyS_Nominal = l.py()*shiftP;
    double pzS_Nominal = l.pz()*shiftP;
    double massS_Nominal = l.mass()*shiftMass;
    double enS_Nominal = TMath::Sqrt(pxS_Nominal*pxS_Nominal + pyS_Nominal*pyS_Nominal + pzS_Nominal*pzS_Nominal + massS_Nominal*massS_Nominal);
    math::XYZTLorentzVectorD p4S_Nominal( pxS_Nominal, pyS_Nominal, pzS_Nominal, enS_Nominal );

    l.setP4( p4S_Nominal );
    
    //Up and Down (+3/-3%) variations
    const float udShift[2] = {1.03, 0.97}; // 0: UP, 1: DOWN
    // ShiftDown = 0.97;
    // if(NominalUpOrDown=="Nominal") Shift = 1.;
    // if(NominalUpOrDown=="Up") Shift = 1.03;
    // if(NominalUpOrDown=="Down") Shift = 0.97;

    float udshiftP[2] = {1., 1.};
    float udshiftMass[2] = {1., 1.};
    bool isTESShifted = false;
    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. ) {

      isTESShifted = true;

      if(l.decayMode()>=1 && l.decayMode()<=10){
        udshiftP[0] = udShift[0]; // up
        udshiftP[1] = udShift[1]; // down
        udshiftMass[0] = udShift[0]; // up
        udshiftMass[1] = udShift[1]; // down
      }

      else if(l.decayMode()==0){
        udshiftP[0] = udShift[0]; // up
        udshiftP[1] = udShift[1]; // down
        udshiftMass[0] = udshiftMass[1] = 1.; // no mass shift for pi0
      }
      // else if(l.decayMode()==10){
      //   shiftP = Shift;
      //   shiftMass = Shift;
      // }
      else isTESShifted = false;
    }

    if (isTESShifted)
    {
      // up shift
      double pxS = l.px()*udshiftP[0];
      double pyS = l.py()*udshiftP[0];
      double pzS = l.pz()*udshiftP[0];
      double massS = l.mass()*udshiftMass[0];
      double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);

      math::XYZTLorentzVectorD p4SUP( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauUp",p4SUP.px());
      l.addUserFloat("py_TauUp",p4SUP.py());
      l.addUserFloat("pz_TauUp",p4SUP.pz());
      l.addUserFloat("e_TauUp",p4SUP.energy());

      // down shift
      pxS = l.px()*udshiftP[1];
      pyS = l.py()*udshiftP[1];
      pzS = l.pz()*udshiftP[1];
      massS = l.mass()*udshiftMass[1];
      enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4SDOWN( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauDown",p4SDOWN.px());
      l.addUserFloat("py_TauDown",p4SDOWN.py());
      l.addUserFloat("pz_TauDown",p4SDOWN.pz());
      l.addUserFloat("e_TauDown",p4SDOWN.energy());
    }
    l.addUserInt("TauUpExists", isTESShifted ? 1 : 0);
    l.addUserInt("TauDownExists", isTESShifted ? 1 : 0);


    edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
    double totalPx(0);
    if (l.decayMode() == 10 /* or l.decayMode() == 1*/) {
      ///////////////////////////////////////////////////////////////////////////////////////////////
      // Get tracks form PFTau daugthers
      std::vector<reco::TransientTrack> transTrk;
      TransientVertex transVtx;
      // const reco::PFCandidateRefVector cands = l.signalChargedHadrCands();
      reco::CandidatePtrVector 	sigCands = l.signalChargedHadrCands();//signalCands();
     
      edm::Handle<reco::BeamSpot> beamSpot;
      iEvent.getByToken(beamSpotTag, beamSpot);
 

      //Get tracks associated wiht pfPV
      reco::TrackCollection pvTracks;
      TLorentzVector aTrack;
      for(size_t i=0; i<cands->size(); ++i){
	if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
	if(!(*cands)[i].bestTrack()) continue;
    
	unsigned int key = (*cands)[i].vertexRef().key();
	int quality = (*cands)[i].pvAssociationQuality();

	if(key!=0 ||
	   (quality!=pat::PackedCandidate::UsedInFitTight
	    && quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

	pvTracks.push_back(*((*cands)[i].bestTrack()));
      }

      std::vector<reco::TransientTrack> transTracks;  
      //   find  tracks belonging to tau decay

      std::vector<double > tracksToBeRemoved; // compare by Pt due to the conflict of comparing const and not const iterators
      double matchingQuality(0);
      for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr) {
	double deR(999.); 
	double checkqual(0);
	reco::Track closestTrack;
	for(auto iter: pvTracks) {
	  if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end()) continue;
	  if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deR){
	    deR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
	    checkqual=deR;
	    closestTrack = iter;
	  }
	}
	matchingQuality+=checkqual;
	tracksToBeRemoved.push_back(closestTrack.pt());
	transTracks.push_back(transTrackBuilder->build(closestTrack));  //cout<<"  closestTrackiter eta  :  "<<   closestTrack.eta() << "   phi   " << closestTrack.phi() << "    pt  "<< closestTrack.pt() <<endl;
      }
        bool fitOk = false;  
      if(transTracks.size() >= 2 ) {
	AdaptiveVertexFitter avf;
	avf.setWeightThreshold(0.001); 
	try {
	  transVtx = avf.vertex(transTracks, *beamSpot);
	  fitOk = true; 
	} catch (...) {
	  fitOk = false; 
	  std::cout<<"Vtx fit failed!"<<std::endl;
	}
      }

      fitOk = fitOk && transVtx.isValid() && fabs(transVtx.position().x())<1 && fabs(transVtx.position().y())<1;
    
      if(fitOk) {
	///NOTE: we take original vertex z position, as this gives the best reults on CP
	///variables. To be understood; probable reason are missing tracks with Pt<0.95GeV
	SVPos.push_back(transVtx.position().x());
	SVPos.push_back(transVtx.position().y());
	SVPos.push_back(transVtx.position().z());
      }
      reco::Vertex secondaryVertex = transVtx;
 
      SVChi2NDofMatchingQual.push_back(secondaryVertex.chi2());
      SVChi2NDofMatchingQual.push_back(secondaryVertex.ndof());
      SVChi2NDofMatchingQual.push_back(matchingQuality);
      TMatrixTSym<double> svcov(3);
      math::Error<3>::type svCov;
      secondaryVertex.fill(svCov);
      for (int i = 0; i <3; i++){
	for (int j = 0; j < 3; j++) {
	  svcov(i, j) = svCov(i, j);
	  svcov(j, i) = svCov(i, j);
	  // cout<<"  svcov  "<<svcov(i,j)<<endl;
	}
      }
      for (int i = 0; i < 3; i++) {
	for (int j = i; j < 3; j++) {
	  SVCov.push_back(svcov(i, j));
	}
      }
   
      for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr) {
	std::vector<double> iiPionP4;
	iiPionP4.push_back((*itr)->p4().e());
	iiPionP4.push_back((*itr)->p4().px());
	iiPionP4.push_back((*itr)->p4().py());
	iiPionP4.push_back((*itr)->p4().pz());
	iPionCharge.push_back((*itr)->charge());
	iPionP4.push_back(iiPionP4);
	totalPx+=(*itr)->p4().px();
      }
    }
  
    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(l);

    int numChargedParticlesSignalCone = l.signalChargedHadrCands().size();
    int numNeutralHadronsSignalCone = l.signalNeutrHadrCands().size();
    int numPhotonsSignalCone = l.signalGammaCands().size();
    int numParticlesSignalCone = l.signalCands().size();
    int numChargedParticlesIsoCone = l.isolationChargedHadrCands().size();
    int numNeutralHadronsIsoCone = l.isolationNeutrHadrCands().size();
    int numPhotonsIsoCone = l.isolationGammaCands().size();
    int numParticlesIsoCone = l.isolationCands().size();
    float leadChargedParticlePt=l.leadCand()->pt();
    float trackRefPt = (l.leadChargedHadrCand().isNonnull() ? l.leadChargedHadrCand()->pt() : 0.);

    
    //Decay mode
    //int decayMode = -1;
    //int A = l.signalPFChargedHadrCands().size();
    //int B = l.signalPFGammaCands().size();
    //if(A==1&&B==0)decayMode=1;
    //else if(A==1&&B>0)decayMode=2;
    //else if (A==3&&B==0)decayMode=3;
    float tauid = (l.isTauIDAvailable(theDiscriminatorTag) ? l.tauID(theDiscriminatorTag) : -999);
    //printf("A, B, tau %d %d %f \n",A,B,tauid);

    //if(decayMode<0&&tauid==0)edm::LogWarning("TauFiller: Unrecognized decay mode");
    /*

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;
    */
    
    float dxy = 999.;
    float dz  = 999.;
    if (vertexs->size()>0) {
      //dxy = l.dxy();
      //const Vertex* vertex = &(vertexs->front());          
      //dz = l.vertex().z() - vertex[0].z();

      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(l.leadChargedHadrCand().get());
      dz=packedLeadTauCand->dz();
      dxy=packedLeadTauCand->dxy();

      //For some reasons, the reference secondaryVertex() is empty EVEN if hasSecondaryVertex is true
      //To be asked to miniAOD people
     // if(l.hasSecondaryVertex()) {
     //   dz  = l.secondaryVertex().get()->z();    

     // 	cout<<"  secondarry "<< dz <<endl;
     // }
    }
  
   
    //--- Embed user variables
    l.addUserInt("isTESShifted",isTESShifted);
    l.addUserFloat("HPSDiscriminator",tauid); 
    l.addUserFloat("decayMode",l.decayMode()); 
    l.addUserFloat("dxy",dxy); 
    l.addUserFloat("dz",dz); 
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso); 
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso); 
    l.addUserFloat("PFPhotonIso",PFPhotonIso); 
    l.addUserFloat("combRelIsoPF",combRelIsoPF); 
    l.addUserInt("numChargedParticlesSignalCone",numChargedParticlesSignalCone);
    l.addUserInt("numNeutralHadronsSignalCone",numNeutralHadronsSignalCone);
    l.addUserInt("numPhotonsSignalCone",numPhotonsSignalCone);
    l.addUserInt("numParticlesSignalCone",numParticlesSignalCone);
    l.addUserInt("numChargedParticlesIsoCone",numChargedParticlesIsoCone);
    l.addUserInt("numNeutralHadronsIsoCone",numNeutralHadronsIsoCone);
    l.addUserInt("numPhotonsIsoCone",numPhotonsIsoCone);
    l.addUserInt("numParticlesIsoCone",numParticlesIsoCone);
    l.addUserFloat("leadChargedParticlePt",leadChargedParticlePt);
    l.addUserFloat("trackRefPt",trackRefPt); 

    // fill all userfloats
    for (unsigned int iuf = 0; iuf < tauFloatDiscrims_.size(); iuf++)
    {
      string ID = tauFloatDiscrims_.at(iuf);
      l.addUserFloat (ID.c_str(), l.isTauIDAvailable(ID.c_str()) ? l.tauID (ID.c_str()) : -999.);
    }

    // fill all userints
    for (unsigned int iui = 0; iui < tauIntDiscrims_.size(); iui++)
    {
      string ID = tauIntDiscrims_.at(iui);
      int ui = -999;
      if (l.isTauIDAvailable(ID.c_str()))
      {
        ui = ( (l.tauID (ID.c_str()) > 0.5) ? 1 : 0);
      }
      l.addUserInt (ID.c_str(), ui);
    }


    //--- MC parent code 
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    float pxHad=0, pyHad=0, pzHad=0, EHad=0; // hadronic gen tau
    int status=99999, id=99999;

    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();
      status =genL->status();
      id =genL->pdgId();

      //cout << "Tau filler: " << i << " [px, id] = " << l.px() << " , " << l.pdgId() << " | (px, py, pz, e) " << px << " " << py << " " << pz << " " << E << " | ID: " << genL->pdgId() << " | status: " << genL->status() << endl;
   
      // build hadronic gen tau (all visible sons)
      for (unsigned int iDau = 0; iDau < genL->numberOfDaughters(); iDau++)
      {
          const Candidate * Dau = genL->daughter( iDau );
          int dauId = Dau->pdgId();
          if (abs(dauId) != 12 && abs(dauId) != 14 && abs(dauId) != 16)
	    {
              pxHad += Dau->p4().Px();
              pyHad += Dau->p4().Py();
              pzHad += Dau->p4().Pz();
              EHad += Dau->p4().E();
          }
      }


      
      //search if it comes from H
      Handle<edm::View<reco::GenParticle> > genHandle;
      iEvent.getByToken(theGenTag, genHandle);
      for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
        int pdgmot = (&(*genHandle)[ipruned])->pdgId();
        if(abs(pdgmot)==25){
          if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],genL)){
            fromH=1;
            break;
          }
        }
      }
    }
    //    for(unsigned int i =0; i< iPionP4.size(); i++){ cout<<"TauFiller  "<<iPionP4.at(i).size() <<  "pions energy   " << iPionP4.at(i).at(0)<<endl;    }
    l.addUserData<std::vector<double > >( "SVPos", SVPos );
    l.addUserData<std::vector<double > >( "SVCov", SVCov);
    l.addUserData<std::vector<double > >( "SVChi2NDofMatchingQual",  SVChi2NDofMatchingQual);
    l.addUserData<std::vector<std::vector<double > > >( "iPionP4", iPionP4 );
    l.addUserData<std::vector<double > >( "iPionCharge",  iPionCharge);
    l.addUserFloat("genPx",px);
    l.addUserFloat("genPy",py);
    l.addUserFloat("genPz",pz);
    l.addUserFloat("genE",E);
    l.addUserInt("status", status);
    l.addUserInt("id", id);

    l.addUserFloat("fromH",fromH);

    l.addUserFloat("genHadPx",px);
    l.addUserFloat("genHadPy",py);
    l.addUserFloat("genHadPz",pz);
    l.addUserFloat("genHadE",E);

    //     MCHistoryTools mch(iEvent);
    //     if (mch.isMC()) {
    //       int MCParentCode = 0;
    //       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
    //       l.addUserFloat("MCParentCode",MCParentCode);
    //     }
    
    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Tau>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    
    result->push_back(l);
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

