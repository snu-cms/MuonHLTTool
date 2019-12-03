// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

//--- for SimHit association
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

class MuonHLTSeedNtupler : public edm::EDAnalyzer
{
public:
  MuonHLTSeedNtupler(const edm::ParameterSet &iConfig);
  virtual ~MuonHLTSeedNtupler() {};

  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
  virtual void endRun(const edm::Run &iRun, const edm::EventSetup &iSetup);

private:
  void Init();
  void Make_Branch();

  void Fill_IterL3TT(const edm::Event &iEvent);
  void Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  // TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken;

  edm::EDGetTokenT< l1t::MuonBxCollection >                  t_L1Muon_;

  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIterL3OISeedsFromL2Muons_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter2IterL3MuonPixelSeeds_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter3IterL3MuonPixelSeeds_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter2IterL3FromL1MuonPixelSeeds_;
  edm::EDGetTokenT< TrajectorySeedCollection >               t_hltIter3IterL3FromL1MuonPixelSeeds_;

  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIterL3OIMuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter0IterL3MuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter2IterL3MuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter3IterL3MuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter0IterL3FromL1MuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter2IterL3FromL1MuonTrack_;
  edm::EDGetTokenT< edm::View<reco::Track> >               t_hltIter3IterL3FromL1MuonTrack_;

  TTree *NThltIterL3OI_;
  TTree *NThltIter0_;
  TTree *NThltIter2_;
  TTree *NThltIter3_;
  TTree *NThltIter0FromL1_;
  TTree *NThltIter2FromL1_;
  TTree *NThltIter3FromL1_;

  class tmpTSOD {
  private:
    uint32_t TSODDetId;
    float TSODPt;
    float TSODX;
    float TSODY;
    float TSODDxdz;
    float TSODDydz;
    float TSODPx;
    float TSODPy;
    float TSODPz;
    float TSODqbp;
    int TSODCharge;
  public:
    void SetTmpTSOD(const PTrajectoryStateOnDet TSODIn) {
      TSODDetId = TSODIn.detId();
      TSODPt = TSODIn.pt();
      TSODX = TSODIn.parameters().position().x();
      TSODY = TSODIn.parameters().position().y();
      TSODDxdz = TSODIn.parameters().dxdz();
      TSODDydz = TSODIn.parameters().dydz();
      TSODPx = TSODIn.parameters().momentum().x();
      TSODPy = TSODIn.parameters().momentum().y();
      TSODPz = TSODIn.parameters().momentum().z();
      TSODqbp = TSODIn.parameters().qbp();
      TSODCharge = TSODIn.parameters().charge();
    }

    tmpTSOD(const PTrajectoryStateOnDet TSODIn) { SetTmpTSOD(TSODIn); }

    bool operator==(const tmpTSOD& other) const {
      return (
        this->TSODDetId == other.TSODDetId &&
        this->TSODPt == other.TSODPt &&
        this->TSODX == other.TSODX &&
        this->TSODY == other.TSODX &&
        this->TSODDxdz == other.TSODDxdz &&
        this->TSODDydz == other.TSODDydz &&
        this->TSODPx == other.TSODPx &&
        this->TSODPy == other.TSODPy &&
        this->TSODPz == other.TSODPz &&
        this->TSODqbp == other.TSODqbp &&
        this->TSODCharge == other.TSODCharge
      );
    }

    bool operator<(const tmpTSOD& other) const {
      return (this->TSODPt!=other.TSODPt) ? (this->TSODPt < other.TSODPt) : (this->TSODDetId < other.TSODDetId);
    }
  };

  class trkTemplate {
  private:
    int nTrks;
    std::vector<double> trkPts;
    std::vector<double> trkEtas;
    std::vector<double> trkPhis;
    std::vector<int> trkCharges;
    std::vector<int> linkToL3s;
    std::vector<float> bestMatchTP_charge;
    std::vector<int> bestMatchTP_pdgId;
    std::vector<double> bestMatchTP_energy;
    std::vector<double> bestMatchTP_pt;
    std::vector<double> bestMatchTP_eta;
    std::vector<double> bestMatchTP_phi;
    std::vector<double> bestMatchTP_parentVx;
    std::vector<double> bestMatchTP_parentVy;
    std::vector<double> bestMatchTP_parentVz;
    std::vector<int> bestMatchTP_status;
    std::vector<int> bestMatchTP_numberOfHits;
    std::vector<int> bestMatchTP_numberOfTrackerHits;
    std::vector<int> bestMatchTP_numberOfTrackerLayers;
    std::vector<double> bestMatchTP_sharedFraction;
    std::vector<int> matchedTPsize;
  public:
    void clear() {
      nTrks = 0;
      trkPts.clear();
      trkEtas.clear();
      trkPhis.clear();
      trkCharges.clear();
      linkToL3s.clear();
      bestMatchTP_charge.clear();
      bestMatchTP_pdgId.clear();
      bestMatchTP_energy.clear();
      bestMatchTP_pt.clear();
      bestMatchTP_eta.clear();
      bestMatchTP_phi.clear();
      bestMatchTP_parentVx.clear();
      bestMatchTP_parentVy.clear();
      bestMatchTP_parentVz.clear();
      bestMatchTP_status.clear();
      bestMatchTP_numberOfHits.clear();
      bestMatchTP_numberOfTrackerHits.clear();
      bestMatchTP_numberOfTrackerLayers.clear();
      bestMatchTP_sharedFraction.clear();
      matchedTPsize.clear();

      return;
    }

    void setBranch(TTree* tmpntpl, TString name) {
      tmpntpl->Branch("n"+name, &nTrks);
      tmpntpl->Branch(name+"_pt", &trkPts);
      tmpntpl->Branch(name+"_eta", &trkEtas);
      tmpntpl->Branch(name+"_phi", &trkPhis);
      tmpntpl->Branch(name+"_charge", &trkCharges);
      tmpntpl->Branch(name+"_matchedL3", &linkToL3s);
      tmpntpl->Branch(name+"_bestMatchTP_charge", &bestMatchTP_charge);
      tmpntpl->Branch(name+"_bestMatchTP_pdgId", &bestMatchTP_pdgId);
      tmpntpl->Branch(name+"_bestMatchTP_energy", &bestMatchTP_energy);
      tmpntpl->Branch(name+"_bestMatchTP_pt", &bestMatchTP_pt);
      tmpntpl->Branch(name+"_bestMatchTP_eta", &bestMatchTP_eta);
      tmpntpl->Branch(name+"_bestMatchTP_phi", &bestMatchTP_phi);
      tmpntpl->Branch(name+"_bestMatchTP_parentVx", &bestMatchTP_parentVx);
      tmpntpl->Branch(name+"_bestMatchTP_parentVy", &bestMatchTP_parentVy);
      tmpntpl->Branch(name+"_bestMatchTP_parentVz", &bestMatchTP_parentVz);
      tmpntpl->Branch(name+"_bestMatchTP_status", &bestMatchTP_status);
      tmpntpl->Branch(name+"_bestMatchTP_numberOfHits", &bestMatchTP_numberOfHits);
      tmpntpl->Branch(name+"_bestMatchTP_numberOfTrackerHits", &bestMatchTP_numberOfTrackerHits);
      tmpntpl->Branch(name+"_bestMatchTP_numberOfTrackerLayers", &bestMatchTP_numberOfTrackerLayers);
      tmpntpl->Branch(name+"_bestMatchTP_sharedFraction", &bestMatchTP_sharedFraction);
      tmpntpl->Branch(name+"_matchedTPsize", &matchedTPsize);

      return;
    }

    void fill(const reco::Track trk) {
      trkPts.push_back(trk.pt());
      trkEtas.push_back(trk.eta());
      trkPhis.push_back(trk.phi());
      trkCharges.push_back(trk.charge());
      nTrks++;

      return;
    }

    void fillBestTP(const TrackingParticleRef TP) {
      bestMatchTP_charge.push_back(TP->charge());
      bestMatchTP_pdgId.push_back(TP->pdgId());
      bestMatchTP_energy.push_back(TP->energy());
      bestMatchTP_pt.push_back(TP->pt());
      bestMatchTP_eta.push_back(TP->eta());
      bestMatchTP_phi.push_back(TP->phi());
      bestMatchTP_parentVx.push_back(TP->vx());
      bestMatchTP_parentVy.push_back(TP->vy());
      bestMatchTP_parentVz.push_back(TP->vz());
      bestMatchTP_status.push_back(TP->status());
      bestMatchTP_numberOfHits.push_back(TP->numberOfHits());
      bestMatchTP_numberOfTrackerHits.push_back(TP->numberOfTrackerHits());
      bestMatchTP_numberOfTrackerLayers.push_back(TP->numberOfTrackerLayers());

      return;
    }

    void fillDummyTP() {
      bestMatchTP_charge.push_back(-99999.);
      bestMatchTP_pdgId.push_back(-99999);
      bestMatchTP_energy.push_back(-99999.);
      bestMatchTP_pt.push_back(-99999.);
      bestMatchTP_eta.push_back(-99999.);
      bestMatchTP_phi.push_back(-99999.);
      bestMatchTP_parentVx.push_back(-99999.);
      bestMatchTP_parentVy.push_back(-99999.);
      bestMatchTP_parentVz.push_back(-99999.);
      bestMatchTP_status.push_back(-99999);
      bestMatchTP_numberOfHits.push_back(-99999);
      bestMatchTP_numberOfTrackerHits.push_back(-99999);
      bestMatchTP_numberOfTrackerLayers.push_back(-99999);

      return;
    }

    void linkIterL3(int linkNo) { linkToL3s.push_back(linkNo); }
    int matchedIDpassedL3(int idx) { return linkToL3s.at(idx); }
    void fillBestTPsharedFrac(double frac) { bestMatchTP_sharedFraction.push_back(frac); }
    void fillmatchedTPsize(int TPsize) { matchedTPsize.push_back(TPsize); }

    int get_bestMatchTP_pdgId(int idx) { return bestMatchTP_pdgId.at(idx); }
    int get_matchedTPsize(int idx) { return matchedTPsize.at(idx); }

    void print() {
      std::cout << "\nnTrks: " << nTrks << std::endl;
      std::cout << "\t trkPts: " << trkPts.size() << std::endl;
      std::cout << "\t trkEtas: " << trkEtas.size() << std::endl;
      std::cout << "\t trkPhis: " << trkPhis.size() << std::endl;
      std::cout << "\t trkCharges: " << trkCharges.size() << std::endl;
      std::cout << "\t linkToL3s: " << linkToL3s.size() << std::endl;
      std::cout << "\t bestMatchTP_charge: " << bestMatchTP_charge.size() << std::endl;
      std::cout << "\t bestMatchTP_pdgId: " << bestMatchTP_pdgId.size() << std::endl;
      std::cout << "\t bestMatchTP_energy: " << bestMatchTP_energy.size() << std::endl;
      std::cout << "\t bestMatchTP_pt: " << bestMatchTP_pt.size() << std::endl;
      std::cout << "\t bestMatchTP_eta: " << bestMatchTP_eta.size() << std::endl;
      std::cout << "\t bestMatchTP_phi: " << bestMatchTP_phi.size() << std::endl;
      std::cout << "\t bestMatchTP_parentVx: " << bestMatchTP_parentVx.size() << std::endl;
      std::cout << "\t bestMatchTP_parentVy: " << bestMatchTP_parentVy.size() << std::endl;
      std::cout << "\t bestMatchTP_parentVz: " << bestMatchTP_parentVz.size() << std::endl;
      std::cout << "\t bestMatchTP_status: " << bestMatchTP_status.size() << std::endl;
      std::cout << "\t bestMatchTP_numberOfHits: " << bestMatchTP_numberOfHits.size() << std::endl;
      std::cout << "\t bestMatchTP_numberOfTrackerHits: " << bestMatchTP_numberOfTrackerHits.size() << std::endl;
      std::cout << "\t bestMatchTP_numberOfTrackerLayers: " << bestMatchTP_numberOfTrackerLayers.size() << std::endl;
      std::cout << "\t bestMatchTP_sharedFraction: " << bestMatchTP_sharedFraction.size() << std::endl;
      std::cout << "\t matchedTPsize: " << matchedTPsize.size() << std::endl;

      return;
    }
  };

  class seedTemplate {
  private:
    int dir_;
    uint32_t tsos_detId_;
    float tsos_pt_;
    float tsos_pt_val_;
    float tsos_eta_;
    float tsos_phi_;
    float tsos_glob_x_;
    float tsos_glob_y_;
    float tsos_glob_z_;
    int tsos_hasErr_;
    float tsos_err0_;
    float tsos_err1_;
    float tsos_err2_;
    float tsos_err3_;
    float tsos_err4_;
    float tsos_err5_;
    float tsos_err6_;
    float tsos_err7_;
    float tsos_err8_;
    float tsos_err9_;
    float tsos_err10_;
    float tsos_err11_;
    float tsos_err12_;
    float tsos_err13_;
    float tsos_err14_;
    float tsos_x_;
    float tsos_y_;
    float tsos_dxdz_;
    float tsos_dydz_;
    float tsos_px_;
    float tsos_py_;
    float tsos_pz_;
    float tsos_qbp_;
    int tsos_charge_;
    float dR_minDRL1SeedP_;
    float dPhi_minDRL1SeedP_;
    float dR_minDPhiL1SeedX_;
    float dPhi_minDPhiL1SeedX_;
    float dR_minDRL1SeedP_AtVtx_;
    float dPhi_minDRL1SeedP_AtVtx_;
    float dR_minDPhiL1SeedX_AtVtx_;
    float dPhi_minDPhiL1SeedX_AtVtx_;
    int bestMatchTP_pdgId_;
    int matchedTPsize_;
  public:
    void clear() {
      dir_ = -99999;
      tsos_detId_ = 0;
      tsos_pt_ = -99999.;
      tsos_pt_val_ = -99999.;
      tsos_eta_ = -99999.;
      tsos_phi_ = -99999.;
      tsos_glob_x_ = -99999.;
      tsos_glob_y_ = -99999.;
      tsos_glob_z_ = -99999.;
      tsos_hasErr_ = -99999;
      tsos_err0_ = -99999.;
      tsos_err1_ = -99999.;
      tsos_err2_ = -99999.;
      tsos_err3_ = -99999.;
      tsos_err4_ = -99999.;
      tsos_err5_ = -99999.;
      tsos_err6_ = -99999.;
      tsos_err7_ = -99999.;
      tsos_err8_ = -99999.;
      tsos_err9_ = -99999.;
      tsos_err10_ = -99999.;
      tsos_err11_ = -99999.;
      tsos_err12_ = -99999.;
      tsos_err13_ = -99999.;
      tsos_err14_ = -99999.;
      tsos_x_ = -99999.;
      tsos_y_ = -99999.;
      tsos_dxdz_ = -99999.;
      tsos_dydz_ = -99999.;
      tsos_px_ = -99999.;
      tsos_py_ = -99999.;
      tsos_pz_ = -99999.;
      tsos_qbp_ = -99999.;
      tsos_charge_ = -99999;
      dR_minDRL1SeedP_ = -99999.;
      dPhi_minDRL1SeedP_ = -99999.;
      dR_minDPhiL1SeedX_ = -99999.;
      dPhi_minDPhiL1SeedX_ = -99999.;
      dR_minDRL1SeedP_AtVtx_ = -99999.;
      dPhi_minDRL1SeedP_AtVtx_ = -99999.;
      dR_minDPhiL1SeedX_AtVtx_ = -99999.;
      dPhi_minDPhiL1SeedX_AtVtx_ = -99999.;
      bestMatchTP_pdgId_ = -99999;
      matchedTPsize_ = -99999;

      return;
    }

    void setBranch(TTree* tmpntpl) {
      tmpntpl->Branch("dir",          &dir_, "dir/I");
      tmpntpl->Branch("tsos_detId",   &tsos_detId_, "tsos_detId/i");
      tmpntpl->Branch("tsos_pt",      &tsos_pt_, "tsos_pt/F");
      tmpntpl->Branch("tsos_pt_val",  &tsos_pt_val_, "tsos_pt_val/F");
      tmpntpl->Branch("tsos_eta",     &tsos_eta_, "tsos_eta/F");
      tmpntpl->Branch("tsos_phi",     &tsos_phi_, "tsos_phi/F");
      tmpntpl->Branch("tsos_glob_x",  &tsos_glob_x_, "tsos_glob_x/F");
      tmpntpl->Branch("tsos_glob_y",  &tsos_glob_y_, "tsos_glob_y/F");
      tmpntpl->Branch("tsos_glob_z",  &tsos_glob_z_, "tsos_glob_z/F");
      tmpntpl->Branch("tsos_hasErr",  &tsos_hasErr_, "tsos_hasErr/I");
      tmpntpl->Branch("tsos_err0",    &tsos_err0_, "tsos_err0/F");
      tmpntpl->Branch("tsos_err1",    &tsos_err1_, "tsos_err1/F");
      tmpntpl->Branch("tsos_err2",    &tsos_err2_, "tsos_err2/F");
      tmpntpl->Branch("tsos_err3",    &tsos_err3_, "tsos_err3/F");
      tmpntpl->Branch("tsos_err4",    &tsos_err4_, "tsos_err4/F");
      tmpntpl->Branch("tsos_err5",    &tsos_err5_, "tsos_err5/F");
      tmpntpl->Branch("tsos_err6",    &tsos_err6_, "tsos_err6/F");
      tmpntpl->Branch("tsos_err7",    &tsos_err7_, "tsos_err7/F");
      tmpntpl->Branch("tsos_err8",    &tsos_err8_, "tsos_err8/F");
      tmpntpl->Branch("tsos_err9",    &tsos_err9_, "tsos_err9/F");
      tmpntpl->Branch("tsos_err10",   &tsos_err10_, "tsos_err10/F");
      tmpntpl->Branch("tsos_err11",   &tsos_err11_, "tsos_err11/F");
      tmpntpl->Branch("tsos_err12",   &tsos_err12_, "tsos_err12/F");
      tmpntpl->Branch("tsos_err13",   &tsos_err13_, "tsos_err13/F");
      tmpntpl->Branch("tsos_err14",   &tsos_err14_, "tsos_err14/F");
      tmpntpl->Branch("tsos_x",       &tsos_x_, "tsos_x/F");
      tmpntpl->Branch("tsos_y",       &tsos_y_, "tsos_y/F");
      tmpntpl->Branch("tsos_dxdz",    &tsos_dxdz_, "tsos_dxdz/F");
      tmpntpl->Branch("tsos_dydz",    &tsos_dydz_, "tsos_dydz/F");
      tmpntpl->Branch("tsos_px",      &tsos_px_, "tsos_px/F");
      tmpntpl->Branch("tsos_py",      &tsos_py_, "tsos_py/F");
      tmpntpl->Branch("tsos_pz",      &tsos_pz_, "tsos_pz/F");
      tmpntpl->Branch("tsos_qbp",     &tsos_qbp_, "tsos_qbp/F");
      tmpntpl->Branch("tsos_charge",  &tsos_charge_, "tsos_charge/I");
      tmpntpl->Branch("dR_minDRL1SeedP",     &dR_minDRL1SeedP_, "dR_minDRL1SeedP/F");
      tmpntpl->Branch("dPhi_minDRL1SeedP",   &dPhi_minDRL1SeedP_, "dPhi_minDRL1SeedP/F");
      tmpntpl->Branch("dR_minDPhiL1SeedX",   &dR_minDPhiL1SeedX_, "dR_minDPhiL1SeedX/F");
      tmpntpl->Branch("dPhi_minDPhiL1SeedX", &dPhi_minDPhiL1SeedX_, "dPhi_minDPhiL1SeedX/F");
      tmpntpl->Branch("dR_minDRL1SeedP_AtVtx",     &dR_minDRL1SeedP_AtVtx_, "dR_minDRL1SeedP_AtVtx/F");
      tmpntpl->Branch("dPhi_minDRL1SeedP_AtVtx",   &dPhi_minDRL1SeedP_AtVtx_, "dPhi_minDRL1SeedP_AtVtx/F");
      tmpntpl->Branch("dR_minDPhiL1SeedX_AtVtx",   &dR_minDPhiL1SeedX_AtVtx_, "dR_minDPhiL1SeedX_AtVtx/F");
      tmpntpl->Branch("dPhi_minDPhiL1SeedX_AtVtx", &dPhi_minDPhiL1SeedX_AtVtx_, "dPhi_minDPhiL1SeedX_AtVtx/F");
      tmpntpl->Branch("bestMatchTP_pdgId", &bestMatchTP_pdgId_, "bestMatchTP_pdgId/I");
      tmpntpl->Branch("matchedTPsize", &matchedTPsize_, "matchedTPsize/I");

      return;
    }

    void fill(TrajectorySeed seed, edm::ESHandle<TrackerGeometry> tracker) {
      GlobalVector p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
      GlobalPoint x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

      dir_ = seed.direction();
      tsos_detId_ = seed.startingState().detId();
      tsos_pt_ = seed.startingState().pt();
      tsos_pt_val_ = p.perp();
      tsos_eta_ = p.eta();
      tsos_phi_ = p.phi();
      tsos_glob_x_ = x.x();
      tsos_glob_y_ = x.y();
      tsos_glob_z_ = x.z();
      tsos_hasErr_ = seed.startingState().hasError();
      tsos_err0_ = seed.startingState().error(0);
      tsos_err1_ = seed.startingState().error(1);
      tsos_err2_ = seed.startingState().error(2);
      tsos_err3_ = seed.startingState().error(3);
      tsos_err4_ = seed.startingState().error(4);
      tsos_err5_ = seed.startingState().error(5);
      tsos_err6_ = seed.startingState().error(6);
      tsos_err7_ = seed.startingState().error(7);
      tsos_err8_ = seed.startingState().error(8);
      tsos_err9_ = seed.startingState().error(9);
      tsos_err10_ = seed.startingState().error(10);
      tsos_err11_ = seed.startingState().error(11);
      tsos_err12_ = seed.startingState().error(12);
      tsos_err13_ = seed.startingState().error(13);
      tsos_err14_ = seed.startingState().error(14);
      tsos_x_ = seed.startingState().parameters().position().x();
      tsos_y_ = seed.startingState().parameters().position().y();
      tsos_dxdz_ = seed.startingState().parameters().dxdz();
      tsos_dydz_ = seed.startingState().parameters().dydz();
      tsos_px_ = seed.startingState().parameters().momentum().x();
      tsos_py_ = seed.startingState().parameters().momentum().y();
      tsos_pz_ = seed.startingState().parameters().momentum().z();
      tsos_qbp_ = seed.startingState().parameters().qbp();
      tsos_charge_ = seed.startingState().parameters().charge();

      return;
    }

    void fill_L1vars( float dR_minDRL1SeedP,         float dPhi_minDRL1SeedP,
                      float dR_minDPhiL1SeedX ,      float dPhi_minDPhiL1SeedX,
                      float dR_minDRL1SeedP_AtVtx,   float dPhi_minDRL1SeedP_AtVtx,
                      float dR_minDPhiL1SeedX_AtVtx, float dPhi_minDPhiL1SeedX_AtVtx ) {
      dR_minDRL1SeedP_           = dR_minDRL1SeedP;
      dPhi_minDRL1SeedP_         = dPhi_minDRL1SeedP;
      dR_minDPhiL1SeedX_         = dR_minDPhiL1SeedX;
      dPhi_minDPhiL1SeedX_       = dPhi_minDPhiL1SeedX;

      dR_minDRL1SeedP_AtVtx_     = dR_minDRL1SeedP_AtVtx;
      dPhi_minDRL1SeedP_AtVtx_   = dPhi_minDRL1SeedP_AtVtx;
      dR_minDPhiL1SeedX_AtVtx_   = dR_minDPhiL1SeedX_AtVtx;
      dPhi_minDPhiL1SeedX_AtVtx_ = dPhi_minDPhiL1SeedX_AtVtx;

      return;
    }

    void fill_TP( trkTemplate* TTtrack, int index ) {
      // no matched track: bestMatchTP_pdgId_ = -99999, matchedTPsize_ = -99999
      // matched track found but no TPfound: bestMatchTP_pdgId_ = -99999, matchedTPsize_ = 0
      // matched track found and TPfound: proper values are filled
      if( index < 0 )
        return;

      bestMatchTP_pdgId_ = TTtrack->get_bestMatchTP_pdgId(index);
      matchedTPsize_ = TTtrack->get_matchedTPsize(index);
    }

    void fill_ntuple( TTree* tmpntpl ) {
      tmpntpl->Fill();
    }
  };


  std::map<tmpTSOD,unsigned int> hltIterL3OIMuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter0IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter2IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter3IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter0IterL3FromL1MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter2IterL3FromL1MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter3IterL3FromL1MuonTrackMap;

  trkTemplate* TThltIterL3OIMuonTrack = new trkTemplate();
  trkTemplate* TThltIter0IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter2IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter3IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter0IterL3FromL1MuonTrack = new trkTemplate();
  trkTemplate* TThltIter2IterL3FromL1MuonTrack = new trkTemplate();
  trkTemplate* TThltIter3IterL3FromL1MuonTrack = new trkTemplate();

  seedTemplate* ST = new seedTemplate();
  // seedTemplate* SThltIterL3OISeedsFromL2Muons = new seedTemplate();
  // seedTemplate* SThltIter0IterL3MuonPixelSeedsFromPixelTracks = new seedTemplate();
  // seedTemplate* SThltIter2IterL3MuonPixelSeeds = new seedTemplate();
  // seedTemplate* SThltIter3IterL3MuonPixelSeeds = new seedTemplate();
  // seedTemplate* SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = new seedTemplate();
  // seedTemplate* SThltIter2IterL3FromL1MuonPixelSeeds = new seedTemplate();
  // seedTemplate* SThltIter3IterL3FromL1MuonPixelSeeds = new seedTemplate();

  void fill_trackTemplate(const edm::Event &iEvent, edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
    edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_, edm::Handle<TrackingParticleCollection>& TPCollection_,
    std::map<tmpTSOD,unsigned int>& trkMap, trkTemplate* TTtrack);

  void fill_seedTemplate(
  const edm::Event &, edm::EDGetTokenT<TrajectorySeedCollection>&,
  edm::ESHandle<TrackerGeometry>&, std::map<tmpTSOD,unsigned int>&, trkTemplate*, TTree* );
};
