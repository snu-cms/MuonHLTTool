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

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

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

//--- DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"

// -- for L1TkMu propagation
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTrackerBuilder.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"



#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

typedef pair<const DetLayer*, TrajectoryStateOnSurface> LayerTSOS;
typedef pair<const DetLayer*, const TrackingRecHit*> LayerHit;

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

  void Fill_Event(const edm::Event &iEvent);
  void Fill_IterL3TT(const edm::Event &iEvent);
  void Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  // TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken;

  edm::EDGetTokenT< reco::VertexCollection >                 t_offlineVertex_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >         t_PUSummaryInfo_;

  edm::EDGetTokenT< l1t::MuonBxCollection >                  t_L1Muon_;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L2Muon_;

  edm::EDGetTokenT<l1t::TkMuonCollection> t_L1TkMuon_;
  edm::EDGetTokenT<l1t::TkPrimaryVertexCollection> t_L1TkPrimaryVertex_;

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

  edm::EDGetTokenT< reco::GenParticleCollection >            t_genParticle_;

  TTree *NTEvent_;
  TTree *NThltIterL3OI_;
  TTree *NThltIter0_;
  TTree *NThltIter2_;
  TTree *NThltIter3_;
  TTree *NThltIter0FromL1_;
  TTree *NThltIter2FromL1_;
  TTree *NThltIter3FromL1_;

  int runNum_;
  int lumiBlockNum_;
  unsigned long long eventNum_;
  int nVertex_;
  int truePU_;

  int nhltIterL3OI_;
  int nhltIter0_;
  int nhltIter2_;
  int nhltIter3_;
  int nhltIter0FromL1_;
  int nhltIter2FromL1_;
  int nhltIter3FromL1_;

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
    int truePU_;
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
    float dR_minDRL2SeedP_;
    float dPhi_minDRL2SeedP_;
    float dR_minDPhiL2SeedX_;
    float dPhi_minDPhiL2SeedX_;
    float dR_L1TkMuSeedP_;
    float dPhi_L1TkMuSeedP_;
    int bestMatchTP_pdgId_;
    int matchedTPsize_;
    float gen_pt_;
    float gen_eta_;
    float gen_phi_;
  public:
    virtual ~seedTemplate() {}

    void clear_base() {
      truePU_ = -99999;
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
      dR_minDRL2SeedP_ = -99999.;
      dPhi_minDRL2SeedP_ = -99999.;
      dR_minDPhiL2SeedX_ = -99999.;
      dPhi_minDPhiL2SeedX_ = -99999.;
      dR_L1TkMuSeedP_ = -99999.;
      dPhi_L1TkMuSeedP_ = -99999.;
      bestMatchTP_pdgId_ = -99999;
      matchedTPsize_ = -99999;
      gen_pt_ = -99999.;
      gen_eta_ = -99999.;
      gen_phi_ = -99999.;

      return;
    }

    virtual void clear() { clear_base(); }

    void setBranch_base(TTree* tmpntpl) {
      tmpntpl->Branch("truePU",       &truePU_, "truePU/I");
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
      tmpntpl->Branch("dR_minDRL2SeedP",     &dR_minDRL2SeedP_, "dR_minDRL2SeedP/F");
      tmpntpl->Branch("dPhi_minDRL2SeedP",   &dPhi_minDRL2SeedP_, "dPhi_minDRL2SeedP/F");
      tmpntpl->Branch("dR_minDPhiL2SeedX",   &dR_minDPhiL2SeedX_, "dR_minDPhiL2SeedX/F");
      tmpntpl->Branch("dPhi_minDPhiL2SeedX", &dPhi_minDPhiL2SeedX_, "dPhi_minDPhiL2SeedX/F");
      tmpntpl->Branch("dR_L1TkMuSeedP",     &dR_L1TkMuSeedP_,   "dR_L1TkMuSeedP/F");
      tmpntpl->Branch("dPhi_L1TkMuSeedP",   &dPhi_L1TkMuSeedP_, "dPhi_L1TkMuSeedP/F");
      tmpntpl->Branch("bestMatchTP_pdgId", &bestMatchTP_pdgId_, "bestMatchTP_pdgId/I");
      tmpntpl->Branch("matchedTPsize", &matchedTPsize_, "matchedTPsize/I");
      tmpntpl->Branch("gen_pt",      &gen_pt_, "gen_pt/F");
      tmpntpl->Branch("gen_eta",     &gen_eta_, "gen_eta/F");
      tmpntpl->Branch("gen_phi",     &gen_phi_, "gen_phi/F");

      return;
    }

    virtual void setBranch(TTree* tmpntpl) { setBranch_base(tmpntpl); }

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

    void fill_PU( int truePU ) {
      truePU_           = truePU;
      return;
    }

    void fill_L1TkMuvars( float dR_L1TkMuSeedP, float dPhi_L1TkMuSeedP ) {
      dR_L1TkMuSeedP_   = dR_L1TkMuSeedP;
      dPhi_L1TkMuSeedP_ = dPhi_L1TkMuSeedP;
      return;
    }

    void fill_Genvars( float gen_pt, float gen_eta, float gen_phi ) {
      gen_pt_  = gen_pt;
      gen_eta_ = gen_eta;
      gen_phi_ = gen_phi;
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

    void fill_L2vars( float dR_minDRL2SeedP,         float dPhi_minDRL2SeedP,
                      float dR_minDPhiL2SeedX ,      float dPhi_minDPhiL2SeedX ) {
      dR_minDRL2SeedP_           = dR_minDRL2SeedP;
      dPhi_minDRL2SeedP_         = dPhi_minDRL2SeedP;
      dR_minDPhiL2SeedX_         = dR_minDPhiL2SeedX;
      dPhi_minDPhiL2SeedX_       = dPhi_minDPhiL2SeedX;

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

  void fill_trackTemplate(const edm::Event &iEvent, edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
    edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_, edm::Handle<TrackingParticleCollection>& TPCollection_,
    std::map<tmpTSOD,unsigned int>& trkMap, trkTemplate* TTtrack);

  void fill_seedTemplate(
  const edm::Event &, edm::EDGetTokenT<TrajectorySeedCollection>&,
  edm::ESHandle<TrackerGeometry>&, std::map<tmpTSOD,unsigned int>&, trkTemplate*, TTree*, int &nSeed, edm::ESHandle<MagneticField> magfieldH, const edm::EventSetup &iSetup, GeometricSearchTracker* geomTracker );

  // HERE

  class seedL1TSOSTemplate : public seedTemplate {
  private:
    float l1x1_;
    float l1y1_;
    float l1z1_;
    float l1x2_;
    float l1y2_;
    float l1z2_;
    float hitx1_;
    float hity1_;
    float hitz1_;
    float hitx2_;
    float hity2_;
    float hitz2_;
    float l1x3_;
    float l1y3_;
    float l1z3_;
    float hitx3_;
    float hity3_;
    float hitz3_;
    float l1x4_;
    float l1y4_;
    float l1z4_;
    float hitx4_;
    float hity4_;
    float hitz4_;
    int nHits_;

  public:
    ~seedL1TSOSTemplate() {}

    void clearL1Hit_12() {
      l1x1_ = -99999.;
      l1y1_ = -99999.;
      l1z1_ = -99999.;
      l1x2_ = -99999.;
      l1y2_ = -99999.;
      l1z2_ = -99999.;
      hitx1_ = -99999.;
      hity1_ = -99999.;
      hitz1_ = -99999.;
      hitx2_ = -99999.;
      hity2_ = -99999.;
      hitz2_ = -99999.;
      nHits_ = -99999;
    }

    void clearL1Hit_3() {
      l1x3_ = -99999.;
      l1y3_ = -99999.;
      l1z3_ = -99999.;
      hitx3_ = -99999.;
      hity3_ = -99999.;
      hitz3_ = -99999.;
    }

    void clearL1Hit_4() {
      l1x4_ = -99999.;
      l1y4_ = -99999.;
      l1z4_ = -99999.;
      hitx4_ = -99999.;
      hity4_ = -99999.;
      hitz4_ = -99999.;
    }

    void clear() {
      clear_base();
      clearL1Hit_12();
      clearL1Hit_3();
      clearL1Hit_4();
    }

    void setBranch_12(TTree* tmpntpl) {
      tmpntpl->Branch("l1x1", &l1x1_, "l1x1/F");
      tmpntpl->Branch("l1y1", &l1y1_, "l1y1/F");
      tmpntpl->Branch("l1z1", &l1z1_, "l1z1/F");
      tmpntpl->Branch("hitx1", &hitx1_, "hitx1/F");
      tmpntpl->Branch("hity1", &hity1_, "hity1/F");
      tmpntpl->Branch("hitz1", &hitz1_, "hitz1/F");
      tmpntpl->Branch("l1x2", &l1x2_, "l1x2F");
      tmpntpl->Branch("l1y2", &l1y2_, "l1y2F");
      tmpntpl->Branch("l1z2", &l1z2_, "l1z2F");
      tmpntpl->Branch("hitx2", &hitx2_, "hitx2/F");
      tmpntpl->Branch("hity2", &hity2_, "hity2/F");
      tmpntpl->Branch("hitz2", &hitz2_, "hitz2/F");
      tmpntpl->Branch("nHits", &nHits_, "nHits/I");
    }

    void setBranch_3(TTree* tmpntpl) {
      tmpntpl->Branch("l1x3", &l1x3_, "l1x3/F");
      tmpntpl->Branch("l1y3", &l1y3_, "l1y3/F");
      tmpntpl->Branch("l1z3", &l1z3_, "l1z3/F");
      tmpntpl->Branch("hitx3", &hitx3_, "hitx3/F");
      tmpntpl->Branch("hity3", &hity3_, "hity3/F");
      tmpntpl->Branch("hitz3", &hitz3_, "hitz3/F");
    }

    void setBranch_4(TTree* tmpntpl) {
      tmpntpl->Branch("l1x4", &l1x4_, "l1x4/F");
      tmpntpl->Branch("l1y4", &l1y4_, "l1y4/F");
      tmpntpl->Branch("l1z4", &l1z4_, "l1z4/F");
      tmpntpl->Branch("hitx4", &hitx4_, "hitx4/F");
      tmpntpl->Branch("hity4", &hity4_, "hity4/F");
      tmpntpl->Branch("hitz4", &hitz4_, "hitz4/F");
    }

    void setBranch(TTree* tmpntpl) {
      setBranch_base(tmpntpl);
      setBranch_12(tmpntpl);
      setBranch_3(tmpntpl);
      setBranch_4(tmpntpl);
    }

    void fill_12(pair<LayerHit, LayerTSOS> firstHit, pair<LayerHit, LayerTSOS> secondHit, int nHits) {
      auto hit1   = firstHit.first.second;
      auto tsos1  = firstHit.second.second;

      l1x1_ = tsos1.globalPosition().x();
      l1y1_ = tsos1.globalPosition().y();
      l1z1_ = tsos1.globalPosition().z();
      hitx1_ = hit1->globalPosition().x();
      hity1_ = hit1->globalPosition().y();
      hitz1_ = hit1->globalPosition().z();

      auto hit2   = secondHit.first.second;
      auto tsos2  = secondHit.second.second;

      l1x2_ = tsos2.globalPosition().x();
      l1y2_ = tsos2.globalPosition().y();
      l1z2_ = tsos2.globalPosition().z();
      hitx2_ = hit2->globalPosition().x();
      hity2_ = hit2->globalPosition().y();
      hitz2_ = hit2->globalPosition().z();

      nHits_ = nHits;
    }

    void fill_3(pair<LayerHit, LayerTSOS> thirdHit) {
      auto hit3   = thirdHit.first.second;
      auto tsos3  = thirdHit.second.second;

      l1x3_ = tsos3.globalPosition().x();
      l1y3_ = tsos3.globalPosition().y();
      l1z3_ = tsos3.globalPosition().z();
      hitx3_ = hit3->globalPosition().x();
      hity3_ = hit3->globalPosition().y();
      hitz3_ = hit3->globalPosition().z();
    }

    void fill_4(pair<LayerHit, LayerTSOS> fourthHit) {
      auto hit4   = fourthHit.first.second;
      auto tsos4  = fourthHit.second.second;

      l1x4_ = tsos4.globalPosition().x();
      l1y4_ = tsos4.globalPosition().y();
      l1z4_ = tsos4.globalPosition().z();
      hitx4_ = hit4->globalPosition().x();
      hity4_ = hit4->globalPosition().y();
      hitz4_ = hit4->globalPosition().z();
    }
  };

  seedL1TSOSTemplate* theSeeds;

  void testRun(
    const edm::Event &, const edm::EventSetup&,
    edm::EDGetTokenT<TrajectorySeedCollection>&
  );

  vector< LayerTSOS > getTsosOnPixels(
    l1t::TkMuon,
    edm::ESHandle<MagneticField>&,
    const Propagator&,
    GeometricSearchTracker*
  );

  vector< pair<LayerHit, LayerTSOS> > getHitTsosPairs(
    TrajectorySeed,
    edm::Handle<l1t::TkMuonCollection>,
    edm::ESHandle<MagneticField>&,
    const Propagator&,
    GeometricSearchTracker*
  );
};
