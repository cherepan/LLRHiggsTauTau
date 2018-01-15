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
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/Particle.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/TrackParticle.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/ParticleBuilder.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/LorentzVectorParticle.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/PDGInfo.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/TrackHelixVertexFitter.h"
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;


struct sv_pair {
	double flightLength;
	double flightLengthSignificance;
};

TRandom3 *r3 = new TRandom3();

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> ROOT_TTree_vector3D;

struct sv_pair geometrical_SV(
	ROOT_TTree_vector3D& b_1, ROOT_TTree_vector3D& tr1,
	ROOT_TTree_vector3D& b_2, ROOT_TTree_vector3D& tr2,
	ROOT_TTree_vector3D& b_3, ROOT_TTree_vector3D& tr3
	)
	{
	//Float_t tracker_error = 0.002; // approximately systematic error on positions
	// it will cancel out with weights

	TVector3 b_vec1, b_vec2, b_vec3;
	//b_vec1.SetXYZ(-b1x, -b1y, -b1z); // 100% known that z here has giant error -- need to do something with it
	//b_vec2.SetXYZ(-b2x, -b2y, -b2z);
	//b_vec3.SetXYZ(-b3x, -b3y, -b3z);
	b_vec1.SetXYZ(b_1.X(), b_1.Y(), b_1.Z());
	b_vec2.SetXYZ(b_2.X(), b_2.Y(), b_2.Z());
	b_vec3.SetXYZ(b_3.X(), b_3.Y(), b_3.Z());

	// I need just the direction of tracks for geometry
	// thus making copy
	TVector3 t1, t2, t3;
	t1.SetXYZ(tr1.X(), tr1.Y(), tr1.Z());
	t2.SetXYZ(tr2.X(), tr2.Y(), tr2.Z());
	t3.SetXYZ(tr3.X(), tr3.Y(), tr3.Z());
	//t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	//t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	//t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	// weighted bis direction -- used in simple b SV (Friday result)
	TVector3 t_sum = t1 + t2 + t3;

	// root throws warning "zero vector can't be streched"
	// crab jobs crash with it
	// protective programming follows
	struct sv_pair sv_zeros = {0., 0.};
	if (t_sum.Mag() == 0)
		return sv_zeros;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	if (t1.Mag() == 0 || t2.Mag() == 0 || t3.Mag() == 0)
		return sv_zeros;
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	//TVector3 tau = t1+t2+t3;
	//TVector3 tau;
	//tau.SetPtEtaPhi(taupt, taueta, tauphi);
	// tests show that tau direction and the sum are practically the same

	// find the "optimal direction"
	// -- direction of minimal angles betwee tracks and b-s in perpendicular plane

	// 2)
	// just shift bis direction randomly in max phi max theta deviations
	// choose best position, i.e. max sum b-track angles in transverse plane
	// thus, no Z changes
	double max_angle_sum = 0;
	TVector3 max_average = t_sum; // initial best direction is bis

	// find max phi and theta dev around bis
	double max_dPhi = 0, max_dTheta = 0;

	double dPhi = abs(t_sum.Phi() - t1.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;
	dPhi = abs(t_sum.Phi() - t2.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;
	dPhi = abs(t_sum.Phi() - t3.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;

	double dTheta = abs(t_sum.Theta() - t1.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;
	dTheta = abs(t_sum.Theta() - t2.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;
	dTheta = abs(t_sum.Theta() - t3.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;

	for (unsigned int i = 0; i<1000; i++)
		{
		//// uniform search around bis dir +- max dphi
		//double dPhi_shift   = max_dPhi * r3->Uniform() * 2 - max_dPhi;
		//double dTheta_shift = max_dTheta * r3->Uniform() * 2 - max_dTheta;
		//TVector3 direction = t_sum;
		//direction.SetPhi(t_sum.Phi() + dPhi_shift);
		//direction.SetTheta(t_sum.Theta() + dTheta_shift);

		// Gaussian + Markov walk from bis dir
		double dPhi_shift   = r3->Gaus(0, max_dPhi);
		double dTheta_shift = r3->Gaus(0, max_dTheta);
		// shift around current best (in principle I should also reduce sigma..)
		TVector3 direction = max_average;
		direction.SetPhi(max_average.Phi() + dPhi_shift);
		direction.SetTheta(max_average.Theta() + dTheta_shift);

		if (direction.Mag() == 0)
			return sv_zeros;
		direction.SetMag(1); // just in case

		// and to the direction
		// find perpendicular b-s
		TVector3 b_long1 = direction * (b_vec1.Dot(direction));
		TVector3 b_perp1 = b_vec1 - b_long1;
		TVector3 b_long2 = direction * (b_vec2.Dot(direction));
		TVector3 b_perp2 = b_vec2 - b_long2;
		TVector3 b_long3 = direction * (b_vec3.Dot(direction));
		TVector3 b_perp3 = b_vec3 - b_long3;

		// perpendicular parts of tracks
		TVector3 t1_long = direction * (t1.Dot(direction));
		TVector3 t1_perp = t1 - t1_long;
		TVector3 t2_long = direction * (t2.Dot(direction));
		TVector3 t2_perp = t2 - t2_long;
		TVector3 t3_long = direction * (t3.Dot(direction));
		TVector3 t3_perp = t3 - t3_long;

		double angle_sum = b_perp1.Angle(t1_perp) + b_perp2.Angle(t2_perp) + b_perp3.Angle(t3_perp);
		if (angle_sum > max_angle_sum)
			{
			max_angle_sum = angle_sum;
			max_average = direction;
			}
		}

	// and to optimal direction
	// find perpendicular b-s
	TVector3 b_long1 = max_average * (b_vec1.Dot(max_average));
	TVector3 b_perp1 = b_vec1 - b_long1;
	TVector3 b_long2 = max_average * (b_vec2.Dot(max_average));
	TVector3 b_perp2 = b_vec2 - b_long2;
	TVector3 b_long3 = max_average * (b_vec3.Dot(max_average));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// perpendicular parts of tracks
	TVector3 t1_long = max_average * (t1.Dot(max_average));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = max_average * (t2.Dot(max_average));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = max_average * (t3.Dot(max_average));
	TVector3 t3_perp = t3 - t3_long;

	// project found b-s to perp tracks
	// in principle it should not be needed, since the direction is found to fit them together well
	// but let's try to get to simple SV best result
	if (t1_perp.Mag() == 0 || t2_perp.Mag() == 0 || t3_perp.Mag() == 0)
		return sv_zeros;
	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	TVector3 b_long_perp1 = t1_perp * (b_perp1.Dot(t1_perp));
	TVector3 b_long_perp2 = t2_perp * (b_perp2.Dot(t2_perp));
	TVector3 b_long_perp3 = t3_perp * (b_perp3.Dot(t3_perp));

	// [let's try without these for now]

	/*
	TVector3 b_long_perp1 = b_perp1;
	TVector3 b_long_perp2 = b_perp2;
	TVector3 b_long_perp3 = b_perp3;
	*/


	// perpendiculars to bis direction, for reference
	// find perpendicular b-s
	TVector3 b_bis_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_bis_perp1 = b_vec1 - b_bis_long1;
	TVector3 b_bis_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_bis_perp2 = b_vec2 - b_bis_long2;
	TVector3 b_bis_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_bis_perp3 = b_vec3 - b_bis_long3;

	// perpendicular parts of tracks
	TVector3 t1_bis_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_bis_perp = t1 - t1_bis_long;
	TVector3 t2_bis_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_bis_perp = t2 - t2_bis_long;
	TVector3 t3_bis_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_bis_perp = t3 - t3_bis_long;

	// in the perp plane find b long to tracks
	// -- nope, no additional correction to b-s

	// the best point calculation
	// with just transverse b-s
	//TVector3 dV = t1 - t2;
	//TVector3 dB = b_perp1 - b_perp2;
	//double x12 = dV.Dot(dB) / dV.Mag2();
	//dV = t2 - t3;
	//dB = b_perp2 - b_perp3;
	//double x23 = dV.Dot(dB) / dV.Mag2();
	//dV = t3 - t1;
	//dB = b_perp3 - b_perp1;
	//double x31 = dV.Dot(dB) / dV.Mag2();

	// the best point calculation with projected b-s
	TVector3 dV = t1 - t2;
	TVector3 dB1 = b_long_perp1 - b_long_perp2;
	double x12 = - dV.Dot(dB1) / dV.Mag2();

	dV = t2 - t3;
	TVector3 dB2 = b_long_perp2 - b_long_perp3;
	double x23 = - dV.Dot(dB2) / dV.Mag2();

	dV = t3 - t1;
	TVector3 dB3 = b_long_perp3 - b_long_perp1;
	double x31 = - dV.Dot(dB3) / dV.Mag2();

	TVector3 bp12 = b_long_perp1 + x12 * t1;
	TVector3 bp21 = b_long_perp2 + x12 * t2;
	TVector3 bp_1 = 0.5*(bp12 + bp21);

	TVector3 bp23 = b_long_perp2 + x23 * t2;
	TVector3 bp32 = b_long_perp3 + x23 * t3;
	TVector3 bp_2 = 0.5*(bp23 + bp32);

	TVector3 bp31 = b_long_perp3 + x31 * t3;
	TVector3 bp13 = b_long_perp1 + x31 * t1;
	TVector3 bp_3 = 0.5*(bp31 + bp13);

	TVector3 bp_average = 0.3333*(bp_1 + bp_2 + bp_3);
	TVector3 bp_dev1 = bp_1 - bp_average;
	TVector3 bp_dev2 = bp_2 - bp_average;
	TVector3 bp_dev3 = bp_3 - bp_average;

	/*
	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;
	// not weighted averages
	double x_average = (x12 + x23 + x31) / 3;
	double x_deviation = (pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2)) * 0.3333;
	//double x_dev_syst = x_deviation + syst;
	// weighted average with tracker errors
	//double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	//double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);
	*/

	double convergence_factor = 1;

	//double triang_a = 0, triang_b = 0;
	//const double tan60 = 1.732, sin60 = 0.866;

	double conv1 = (bp12 - bp21).Mag();
	double conv2 = (bp23 - bp32).Mag();
	double conv3 = (bp31 - bp13).Mag();
	double conv_frac1 = conv1/dB1.Mag();
	double conv_frac2 = conv2/dB2.Mag();
	double conv_frac3 = conv3/dB3.Mag();
	double conv_frac_sum = conv_frac1 + conv_frac2 + conv_frac3;
	//	double conv_frac_averaged = sqrt(pow(conv_frac1, 2) + pow(conv_frac2, 2) + pow(conv_frac3, 2));

	// SV out of all penalties
	double flightLength = 0, flightLengthSignificance = 0;
	// by relative convergence volume
	convergence_factor *= 1 / (1 + conv_frac1 * conv_frac2 * conv_frac3 / 0.027); // 0.027 = 0.3*0.3*0.3 -- when fractions are equal
	// by fraction sum-s, extracting correlation of divergences
	convergence_factor *= 1 / (1 + (conv_frac1/conv_frac_sum) * (conv_frac2/conv_frac_sum) * (conv_frac3/conv_frac_sum));

	// sign of flight
	if (bp_average.Dot(t_sum) > 0)
		{
		flightLength = bp_average.Mag();
		flightLengthSignificance = flightLength * convergence_factor / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		}
	else
		{
		flightLength = - bp_average.Mag();
		flightLengthSignificance = flightLength * convergence_factor / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		}

	struct sv_pair SV = {.flightLength = flightLength, .flightLengthSignificance = flightLengthSignificance};
	return SV;
	}







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
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > tracks_;
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
  tracks_(consumes<edm::View<pat::PackedCandidate>> (edm::InputTag("packedPFCandidates"))),
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
  TLorentzVector aTrack;
  reco::TrackCollection pvTracks;
 
     
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


  // TRACKS
  //
  edm::Handle<edm::View<pat::PackedCandidate> > tracksHandle;
  iEvent.getByToken(tracks_, tracksHandle);
  //if (tracksHandle.isValid()) tracks = *tracksHandle;
  const edm::View<pat::PackedCandidate>* track_cands = tracksHandle.product();


  reco::TrackCollection pvertexTracks;
  reco::TrackCollection allTracks; // for taus (with possible SV) (testing now)

  for(size_t i=0; i<track_cands->size(); ++i)
    {
      if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue;
      if(!(*track_cands)[i].bestTrack()) continue;
      
      unsigned int key = (*track_cands)[i].vertexRef().key();
      int quality = (*track_cands)[i].pvAssociationQuality();
      
      // here I need to select "good" tracks
      // save them to all tracks
      // and if they belong to PV save them to pv tracks
      if (!(key!=0 ||
	    (quality!=pat::PackedCandidate::UsedInFitTight
	     && quality!=pat::PackedCandidate::UsedInFitLoose)))// continue;
	{
	  pvertexTracks.push_back(*((*track_cands)[i].bestTrack()));
	  // allTracks.push_back(*((*track_cands)[i].bestTrack())); // test for HelixLine Momentum is zero
	}
      
      // TODO: add requirement of "goodness"?
      allTracks.push_back(*((*track_cands)[i].bestTrack()));
    }



  // Output collection
  auto_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int itau = 0; itau < tauHandle->size(); ++itau){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[itau].get()));
    std::vector<double > SVPos;     
    std::vector<double > SVCov;    
    std::vector<double > PFTau_TIP_PVPos;     
    std::vector<double > PFTau_TIP_PVCov;    


    std::vector<double > PFTauTrackLV;    
    float PFTauTrack_deltaR=999.;
    std::vector<std::vector<double> > iPionP4;
    std::vector<std::vector<double> > iRefitPionP4;
    std::vector<double> iPionCharge;
    std::vector<double> iRefitPionCharge;

    std::vector<double> SVChi2NDofMatchingQual;
    int TauTrackFiller_trackCharge=-999;
    int TauTrackFiller_pdgid=-999;
    float TauTrackFiller_B=-999;
    float TauTrackFiller_M=-999;
    std::vector<double>  TauTrackFiller_par;
    std::vector<double>  TauTrackFiller_cov;
    float GEOMFlightLenght(-999);
    float GEOMFlightLenghtSignificance(-999);


    int a1_charge=-999;
    int a1_pdgid=-999;
    float a1_B=-999.;
    float a1_M=-999.; 
    std::vector<double>  PFTau_a1_lvp;
    std::vector<double>  PFTau_a1_cov;


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
      // else if(l.decayMode()==110){
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

  

  

     // if (l.decayMode() == 0  or l.decayMode() == 1 ) 
     //      {

     // 	   reco::Vertex fakeVertex = vertexs->front();
     // 	   reco::CandidatePtrVector 	sigCands = l.signalChargedHadrCands();
     // 	   //Get tracks associated wiht pfPV
     // 	   reco::TrackCollection pvTracks;
     // 	   TLorentzVector aTrack;
     // 	   for(size_t i=0; i<cands->size(); ++i){
     // 	     if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
     // 	     if(!(*cands)[i].bestTrack()) continue;
    	     
     // 	     unsigned int key = (*cands)[i].vertexRef().key();
     // 	     int quality = (*cands)[i].pvAssociationQuality();
    	     
     // 	     if(key!=0 ||
     // 		(quality!=pat::PackedCandidate::UsedInFitTight
     // 		 && quality!=pat::PackedCandidate::UsedInFitLoose)) continue;
    	     
     // 	     pvTracks.push_back(*((*cands)[i].bestTrack()));
     // 	   }
    	   
     // 	   //---------- find PV track belonging to sigCands
     // 	   double deR(999.); 
     // 	   reco::Track RefToTauTrack;
     // 	   for(auto iter: pvTracks) {
     // 	     //  if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end()) continue;
     // 	     if( sqrt(pow(iter.eta() - l.leadChargedHadrCand()->p4().eta(),2) + pow(iter.phi() - l.leadChargedHadrCand()->p4().phi(),2))  < deR){
     // 	       deR = sqrt(pow(iter.eta() - l.leadChargedHadrCand()->p4().eta(),2) + pow(iter.phi() - l.leadChargedHadrCand()->p4().phi(),2));
     // 	       RefToTauTrack = iter;
     // 	     }
     // 	   }
     // 	   PFTauTrackLV.push_back(l.leadChargedHadrCand()->p4().e());    
     // 	   PFTauTrackLV.push_back(l.leadChargedHadrCand()->p4().px());    
     // 	   PFTauTrackLV.push_back(l.leadChargedHadrCand()->p4().py());    
     // 	   PFTauTrackLV.push_back(l.leadChargedHadrCand()->p4().pz());    

     // 	   //const reco::Track *TauTrack  = (*itr)->bestTrack();
     // 	   const reco::Track *TauTrack  =  &RefToTauTrack;
    	   
     // 	   GlobalPoint pvpoint(TauTrack->vx(), TauTrack->vy(), TauTrack->vz());
     // 	   reco::TransientTrack transTrk = transTrackBuilder->build(TauTrack);
     // 	   TrackParticle tautrackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
     // 	   PFTauTrack_deltaR=deR;

     // 	   if(deR< 0.01){

     // 	     TauTrackFiller_trackCharge=tautrackparticle.Charge();
     // 	     TauTrackFiller_pdgid=tautrackparticle.PDGID();
     // 	     TauTrackFiller_B=tautrackparticle.BField();
     // 	     TauTrackFiller_M=tautrackparticle.Mass();
    	     
     // 	     for (int i = 0; i < tautrackparticle.NParameters(); i++) {
     // 	       TauTrackFiller_par.push_back(tautrackparticle.Parameter(i));
     // 	       for (int j = i; j <tautrackparticle.NParameters(); j++) {
     // 		 TauTrackFiller_cov.push_back(tautrackparticle.Covariance(i, j));
     // 	       }
     // 	     }
     // 	   }else
     // 	         {
     // 		   TauTrackFiller_trackCharge=-999;
     // 		   TauTrackFiller_pdgid=-999;
     // 		   TauTrackFiller_B=-999;
     // 		   TauTrackFiller_M=-999;
     // 		 }
     // 	  }else
     //            {
     // 		 TauTrackFiller_trackCharge=-999;
     // 		 TauTrackFiller_pdgid=-999;
     // 		 TauTrackFiller_B=-999;
     // 		 TauTrackFiller_M=-999;
     // 		}
  
   





    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotTag, beamSpot);
 


    std::vector<reco::TransientTrack> transTracksForPVRefit;  
   
    //   find  tracks belonging to tau decay

    reco::CandidatePtrVector tauproducts = l.signalChargedHadrCands();//signalCands();
    reco::TrackCollection RemovedTracks;
    //    double matchingQuality(0);
    //     std::cout<<"TauFiller:  signalChargedHadrCands: "<<std::endl;
    for (reco::CandidatePtrVector::const_iterator itr = tauproducts.begin(); itr != tauproducts.end(); ++itr) {
      double deR(999.); 
      reco::Track closestTrack;
      //	std::cout<<" px   "<< (*itr)->p4().Px() << " charge   "<<(*itr)->charge() <<std::endl;
      for(auto iter: allTracks) {
	//	if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end()) continue;
	if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deR){
	  deR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
	  closestTrack = iter;
	}
      }
      
      RemovedTracks.push_back(closestTrack);
    }
    
     for(auto iterpv: pvTracks) {
       bool matched(false);
       for(auto iterem: RemovedTracks) {
	 if(sqrt(pow(iterem.eta() - iterpv.eta(),2) + pow(iterem.phi() - iterpv.phi(),2))  < 0.0001) matched =true;
       }
       if(matched)continue;
       transTracksForPVRefit.push_back(transTrackBuilder->build(iterpv));
     }
     bool fitOk = false;  
     TransientVertex transRefitPVVtx;
     if(transTracksForPVRefit.size() >= 2 ) {
       AdaptiveVertexFitter avf;
       avf.setWeightThreshold(0.1); 
       try {
	 transRefitPVVtx = avf.vertex(transTracksForPVRefit, *beamSpot);
	 fitOk = true; 
       } catch (...) {
	 fitOk = false; 
	 std::cout<<"Vtx fit failed!"<<std::endl;
       }
     }



     // if(fitOk){std::cout<<"bla"<<std::endl;}
     if(fitOk ){
       PFTau_TIP_PVPos.push_back(transRefitPVVtx.position().x());
       PFTau_TIP_PVPos.push_back(transRefitPVVtx.position().y());
       PFTau_TIP_PVPos.push_back(transRefitPVVtx.position().z());
       
       reco::Vertex PrimaryTIPVertex = transRefitPVVtx;
       
       TMatrixTSym<double> pvTIPcov(3);
       math::Error<3>::type pvTIPCov;
       PrimaryTIPVertex.fill(pvTIPCov);
       for (int i = 0; i <3; i++){
	 for (int j = 0; j < 3; j++) {
	   pvTIPcov(i, j) = pvTIPCov(i, j);
	   pvTIPcov(j, i) = pvTIPCov(i, j);
	   // cout<<"  svcov  "<<svcov(i,j)<<endl;
	 }
       }
       for (int i = 0; i < 3; i++) {
	 for (int j = i; j < 3; j++) {
	   PFTau_TIP_PVCov.push_back(pvTIPcov(i, j));
	 }
       }
     }

     

     //  std::cout<<"  decay mode    "<<l.decayMode() <<" size  "<< pvTracks.size() << "after removal " <<transTracksForPVRefit.size() <<  "   n tracks to be removed  "<< RemovedTracks.size()<<std::endl;
     //    std::cout<<"  tau sig  "<<l.flightLengthSig()<<" decay Mode  " <<l.decayMode() <<std::endl;
  
    if (l.decayMode() == 10 /* or l.decayMode() == 1*/) {
      ///////////////////////////////////////////////////////////////////////////////////////////////
      // Get tracks form PFTau daugthers
      std::vector<reco::TransientTrack> transTrk;
      TransientVertex transVtx;
      // const reco::PFCandidateRefVector cands = l.signalChargedHadrCands();
      reco::CandidatePtrVector 	sigCands = l.signalChargedHadrCands();//signalCands();
     
      edm::Handle<reco::BeamSpot> beamSpot;
      iEvent.getByToken(beamSpotTag, beamSpot);
 


      std::vector<reco::TransientTrack> transTracks;  
      //   find  tracks belonging to tau decay

      std::vector<double > tracksToBeRemoved; // compare by Pt due to the conflict of comparing const and not const iterators
      double matchingQuality(0);
      //     std::cout<<"TauFiller:  signalChargedHadrCands: "<<std::endl;
      for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr) {
	double deR(999.); 
	double checkqual(0);
	reco::Track closestTrack;
	//	std::cout<<" px   "<< (*itr)->p4().Px() << " charge   "<<(*itr)->charge() <<std::endl;
	for(auto iter: allTracks) {
	  if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end()) continue;
	  if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deR){
	    deR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
	    checkqual=deR;
	    closestTrack = iter;
	  }
	}

	matchingQuality+=checkqual;
	tracksToBeRemoved.push_back(closestTrack.pt());
	if(closestTrack.pt()!=0)transTracks.push_back(transTrackBuilder->build(closestTrack));  //cout<<"  closestTrackiter eta  :  "<<   closestTrack.eta() << "   phi   " << closestTrack.phi() << "    pt  "<< closestTrack.pt() <<endl;
	// std::cout<<" tracks to be reffited  "<<  matchingQuality  <<std::endl;
	//  std::cout<<" px,py,pz,pt  "<< closestTrack.px() << " "<< closestTrack.py() << " "<< closestTrack.pz() << " "<< closestTrack.py() <<std::endl;
      }
        bool fitOk = false;  
	//std::cout<<"  ----  "<<std::endl;

	//pvTracksForRefit

	if(transTracks.size() >= 2 ) {
	  KalmanVertexFitter kvf;
	  try {
	    transVtx = kvf.vertex(transTracks);
	    fitOk = true; 
	  } catch (...) {
	    fitOk = false; 
	    std::cout<<"Vtx fit failed!"<<std::endl;
	  }
	}

      fitOk = fitOk && transVtx.isValid() && fabs(transVtx.position().x())<1 && fabs(transVtx.position().y())<1;
      //      std::cout<<" is fit OK ??  "<<   fitOk <<std::endl;
    
      if(fitOk) {
	///NOTE: we take original vertex z position, as this gives the best reults on CP
	///variables. To be understood; probable reason are missing tracks with Pt<0.95GeV
	SVPos.push_back(transVtx.position().x());
	SVPos.push_back(transVtx.position().y());
	SVPos.push_back(transVtx.position().z());
	//	std::cout<<" SV  "<< transVtx.position().x() <<" " <<transVtx.position().y() <<"   "<< transVtx.position().z() <<std::endl;


	 //--------------------------  geometry calculation of PVSVSignificance. The code was originally written by  A. Toldaiev -------------------------
	 std::vector<Float_t> NT_tau_SV_fit_track_OS_matched_track_dR;
	 std::vector<Float_t> NT_tau_SV_fit_track_SS1_matched_track_dR;
	 std::vector<Float_t> NT_tau_SV_fit_track_SS2_matched_track_dR;
	 std::vector<Float_t> NT_tau_SV_fit_track_matchQ;

	 std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_tau_SV_fit_track_OS_p4;
	 std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_tau_SV_fit_track_SS1_p4;
	 std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_tau_SV_fit_track_SS2_p4;

	  std::vector<Int_t> NT_tau_SV_fit_track_OS_matched_track_vtxkey;
	  std::vector<Int_t> NT_tau_SV_fit_track_SS1_matched_track_vtxkey;
	  std::vector<Int_t> NT_tau_SV_fit_track_SS2_matched_track_vtxkey;
	  std::vector<Int_t> NT_tau_SV_fit_track_OS_matched_track_vtxQ;
	  std::vector<Int_t> NT_tau_SV_fit_track_SS1_matched_track_vtxQ;
	  std::vector<Int_t> NT_tau_SV_fit_track_SS2_matched_track_vtxQ;

	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_OS_matched_track_b;
	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_SS1_matched_track_b;
	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_SS2_matched_track_b;
	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_OS_matched_track_p3;
	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_SS1_matched_track_p3;
	  std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag> >  NT_tau_SV_fit_track_SS2_matched_track_p3;




	
	  //	 std::cout<<" sigCands  "<< sigCands.size()<< "  " <<l.decayMode()<< std::endl;
	 if(sigCands.size() <3)continue;
	 for (reco::CandidatePtrVector::const_iterator itr_cand = sigCands.begin(); itr_cand != sigCands.end(); ++itr_cand)
	   {

	     NT_tau_SV_fit_track_OS_p4.push_back((*itr_cand)->p4());

	   }
	
	 // loop through tracks and save their impact parameters
	 // and match quality dR
	 double min_dR_os(99999.), min_dR_ss1(99999.), min_dR_ss2(99999.);
	 int matched_track_OS = -1, matched_track_SS1 = -1, matched_track_SS2 = -1;
	 for(size_t i=0; i<track_cands->size(); ++i)
	   {
	     // TODO: these requirements are probably the reasone some tracks are not found for tau sigCands
	     if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue;
	     if(!(*track_cands)[i].bestTrack()) continue;
	  
	     auto track = (*track_cands)[i].bestTrack();
	     // if(!track->isNonNull()) continue;
	     // find closest matches to general track	
	     // std::cout<<" deb 10"<<std::endl;     
	     // std::cout<<"NT_tau_SV_fit_track_OS_p4.size()"<<NT_tau_SV_fit_track_OS_p4.size()<<std::endl;     
	     // std::cout<<"NT_tau_SV_fit_track_SS1_p4.size()"<<NT_tau_SV_fit_track_SS1_p4.size()<<std::endl;     
	     // std::cout<<"NT_tau_SV_fit_track_SS2_p4.size()"<<NT_tau_SV_fit_track_SS2_p4.size()<<std::endl;     
 

	     double dR_os  = sqrt(pow(track->eta() - sigCands[0]->p4().eta() , 2) + pow(track->phi() - sigCands[0]->p4().phi(), 2));
	     double dR_ss1 = sqrt(pow(track->eta() - sigCands[1]->p4().eta(), 2) + pow(track->phi() - sigCands[1]->p4().phi(), 2));
	     double dR_ss2 = sqrt(pow(track->eta() - sigCands[2]->p4().eta(), 2) + pow(track->phi() - sigCands[2]->p4().phi(), 2));
	     if (dR_os < min_dR_os)
	       {
	 	min_dR_os = dR_os;
	 	matched_track_OS = i;
	       }
	     if (dR_ss1 < min_dR_ss1)
	       {
	 	min_dR_ss1 = dR_ss1;
	 	matched_track_SS1 = i;
	       }
	     
	     if (dR_ss2 < min_dR_ss2)
	       {
	 	min_dR_ss2 = dR_ss2;
	 	matched_track_SS2 = i;
	       }
	   }
	 // tracks are matched, save parameters
	 //	  std::cout<<"  deb 30 "<<std::endl;
	 // quality of match

	
	 // OS
	 int track_index;
	 track_index = matched_track_OS;

	 auto ref_vertex = *((*track_cands)[track_index].vertexRef());
	 auto closest_point = (*track_cands)[track_index].vertex();
	 auto distance1 = closest_point - ref_vertex.position();
	 // distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
	 //impact.SetXYZ(distance.x(), distance.y(), distance.z());
	
	
	 // SS1
	 track_index = matched_track_SS1;
	
	 ref_vertex = *((*track_cands)[track_index].vertexRef());
	 closest_point = (*track_cands)[track_index].vertex();
	 auto distance2 = closest_point - ref_vertex.position();
	 // distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
	 //impact.SetXYZ(distance.x(), distance.y(), distance.z());
	
	 // SS2
	 track_index = matched_track_SS2;
	
	
	 ref_vertex = *((*track_cands)[track_index].vertexRef());
	 closest_point = (*track_cands)[track_index].vertex();
	 auto distance3 = closest_point - ref_vertex.position();
	 // distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
	 //impact.SetXYZ(distance.x(), distance.y(), distance.z());
	
	
	 //TVector3 tr_ss2;
	 //TVector3 tr_os ;
	 //TVector3 tr_ss1;
	 //tr_ss2.SetXYZ(NT_tau_SV_fit_track_SS2_p4.back().X(), NT_tau_SV_fit_track_SS2_p4.back().Y(), NT_tau_SV_fit_track_SS2_p4.back().Z());
	 //tr_os .SetXYZ(NT_tau_SV_fit_track_OS_p4 .back().X(), NT_tau_SV_fit_track_OS_p4 .back().Y(), NT_tau_SV_fit_track_OS_p4 .back().Z());
	 //tr_ss1.SetXYZ(NT_tau_SV_fit_track_SS1_p4.back().X(), NT_tau_SV_fit_track_SS1_p4.back().Y(), NT_tau_SV_fit_track_SS1_p4.back().Z());
	
	 // let's save also Vector3 of tracks? to not have to convert everything every time
	 NT_tau_SV_fit_track_OS_matched_track_p3.push_back(sigCands[0]->p4().Vect());
	 NT_tau_SV_fit_track_OS_matched_track_p3.push_back(sigCands[1]->p4().Vect());
	 NT_tau_SV_fit_track_OS_matched_track_p3.push_back(sigCands[2]->p4().Vect());
	 // geometrical SV
	 struct sv_pair geom_SV = geometrical_SV(
						 distance1 , NT_tau_SV_fit_track_OS_matched_track_p3.at(0),
					         distance2 , NT_tau_SV_fit_track_OS_matched_track_p3.at(1),
					         distance3 , NT_tau_SV_fit_track_OS_matched_track_p3.at(2));








	 // std::cout<<" deb 60"<<std::endl;
	 GEOMFlightLenght=geom_SV.flightLength;
	 GEOMFlightLenghtSignificance=geom_SV.flightLengthSignificance;

	 //	 std::cout<<" geom_SV.flightLengthSignificance  " <<geom_SV.flightLengthSignificance <<std::endl;
      

    
	 //--------------------------  geometry calculation of PVSVSignificance from A. Toldaiev -------------------------



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

      }

      //      std::cout<<" refiter.size     "<< secondaryVertex.refittedTracks().size()<<std::endl;
      ////////////////////////////////////////////////////////////////////////////////
       if(fitOk &&  transTracks.size() >2) {
       	LorentzVectorParticle a1;
	GlobalPoint sv(secondaryVertex.position().x(), secondaryVertex.position().y(), secondaryVertex.position().z());
	KinematicParticleFactoryFromTransientTrack kinFactory;
       	float piMassSigma(sqrt(pow(10., -12.))), piChi(0.0), piNdf(0.0);
       	std::vector<RefCountedKinematicParticle> pions;
       	for (unsigned int i = 0; i <transTracks.size(); i++)
       	  pions.push_back(kinFactory.particle(transTracks.at(i), PDGInfo::pi_mass(), piChi, piNdf, sv, piMassSigma));
      
       	KinematicParticleVertexFitter kpvFitter;
       	RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
	if(jpTree->isValid()){
	  jpTree->movePointerToTheTop();
	  const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
	  AlgebraicSymMatrix77 cov = jpTree->currentParticle()->currentState().kinematicParametersError().matrix();
	  // get pions
	  double c(0);
	  std::vector<reco::Track> Tracks;
	  std::vector<LorentzVectorParticle> ReFitPions;
	  for (unsigned int i = 0; i < transTracks.size(); i++) {
	    std::vector<double> iPionP4;
	    std::vector<double> iPionCharge;
	    c += transTracks.at(i).charge();
	    ReFitPions.push_back(ParticleBuilder::CreateLorentzVectorParticle(transTracks.at(i), transTrackBuilder, secondaryVertex, true, true));
	    iPionP4.push_back(ReFitPions.at(i).LV().E());
	    iPionP4.push_back(ReFitPions.at(i).LV().Px());
	    iPionP4.push_back(ReFitPions.at(i).LV().Py());
	    iPionP4.push_back(ReFitPions.at(i).LV().Pz());
	    iRefitPionP4.push_back(iPionP4);
	    iRefitPionCharge.push_back(transTracks.at(i).charge());
	    //	    ReFitPions.at(i).LVCov().Print();
	  }
	
	  // now covert a1 into LorentzVectorParticle
	  TMatrixT<double> a1_par(LorentzVectorParticle::NLorentzandVertexPar, 1);
	  TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
	  for (int i = 0; i < LorentzVectorParticle::NLorentzandVertexPar; i++) {
	    a1_par(i, 0) = parameters(i);
	    for (int j = 0; j < LorentzVectorParticle::NLorentzandVertexPar; j++) {
	      a1_cov(i, j) = cov(i, j);
	    }
	  }
	  a1 = LorentzVectorParticle(a1_par, a1_cov, abs(PDGInfo::a_1_plus) * c, c, transTrackBuilder->field()->inInverseGeV(sv).z());
	  a1_charge=a1.Charge();
	  a1_pdgid=a1.PDGID();
	  a1_B=a1.BField();
	  a1_M=a1.Mass();
	  for (int i = 0; i < a1.NParameters(); i++) {
	    PFTau_a1_lvp.push_back(a1.Parameter(i));
	    for (int j = i; j < a1.NParameters(); j++) {
	      PFTau_a1_cov.push_back(a1.Covariance(i, j));
	    }
	  }
	}
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
    l.addUserFloat("FLSig",l.flightLengthSig()); 
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

   l.addUserFloat("GEOMFlightLenght",GEOMFlightLenght); 
    l.addUserFloat("GEOMFlightLenghtSignificance",GEOMFlightLenghtSignificance); 



    //    std::cout<< " Tau " <<itau << std::endl;
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
    l.addUserData<std::vector<std::vector<double > > >( "iRefitPionP4", iRefitPionP4 );
    l.addUserData<std::vector<double > >( "iRefitPionCharge",  iRefitPionCharge);
    l.addUserData<std::vector<double > >( "PFTau_TIP_PVPos", PFTau_TIP_PVPos);
    l.addUserData<std::vector<double > >( "PFTau_TIP_PVCov", PFTau_TIP_PVCov);

    // l.addUserData<std::vector<double > >( "SVPos", SVPos );
    // l.addUserData<std::vector<double > >( "SVCov", SVCov);
  
    l.addUserInt("a1_charge",  a1_charge);
    l.addUserInt("a1_pdgid",  a1_pdgid);
    l.addUserFloat("a1_B",  a1_B);
    l.addUserFloat("a1_M",  a1_M);
    l.addUserData<std::vector<double > >( "PFTau_a1_lvp", PFTau_a1_lvp );
    l.addUserData<std::vector<double > >( "PFTau_a1_cov", PFTau_a1_cov );



    l.addUserInt( "TauTrackFiller_trackCharge", TauTrackFiller_trackCharge);
    l.addUserInt( "TauTrackFiller_pdgid",  TauTrackFiller_pdgid);
    l.addUserFloat( "TauTrackFiller_B",  TauTrackFiller_B);
    l.addUserFloat( "TauTrackFiller_M", TauTrackFiller_M );
    l.addUserData<std::vector<double  > >( "TauTrackFiller_par",  TauTrackFiller_par);
    l.addUserData<std::vector<double  > >( "TauTrackFiller_cov",  TauTrackFiller_cov);

    l.addUserData<std::vector<double  > >( "PFTauTrackLV", PFTauTrackLV );
    l.addUserFloat( "PFTauTrack_deltaR", PFTauTrack_deltaR );

 
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

