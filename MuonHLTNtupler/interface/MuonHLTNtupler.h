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
#include "DataFormats/Common/interface/Ref.h"
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

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
////////////////////////////
// DETECTOR GEOMETRY HEADERS
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

#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

class MuonHLTNtupler : public edm::EDAnalyzer
{
public:
  MuonHLTNtupler(const edm::ParameterSet &iConfig);
  virtual ~MuonHLTNtupler() {};

  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
  virtual void endRun(const edm::Run &iRun, const edm::EventSetup &iSetup);

private:
  void Init();
  void Make_Branch();
  void Fill_L1Track(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void Fill_HLT(const edm::Event &iEvent, bool isMYHLT);
  void Fill_Muon(const edm::Event &iEvent);
  void Fill_HLTMuon(const edm::Event &iEvent);
  void Fill_L1Muon(const edm::Event &iEvent);
  void Fill_GenParticle(const edm::Event &iEvent);

  //For Rerun (Fill_IterL3*)
  void Fill_IterL3(const edm::Event &iEvent);
  void Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  bool SavedTriggerCondition( std::string& pathName );
  bool SavedFilterCondition( std::string& filterName );

  bool isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx);

  bool doSeed;
  bool DebugMode;
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;

  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken;

  edm::EDGetTokenT< std::vector<reco::Muon> >                t_offlineMuon_;
  edm::EDGetTokenT< reco::VertexCollection >                 t_offlineVertex_;
  edm::EDGetTokenT< edm::TriggerResults >                    t_triggerResults_;
  edm::EDGetTokenT< trigger::TriggerEvent >                  t_triggerEvent_;
  edm::EDGetTokenT< edm::TriggerResults >                    t_myTriggerResults_;
  edm::EDGetTokenT< trigger::TriggerEvent >                  t_myTriggerEvent_;

  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L3Muon_;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L2Muon_;
  edm::EDGetTokenT< l1t::MuonBxCollection >                  t_L1Muon_;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_TkMuon_;

  edm::EDGetTokenT< std::vector<reco::MuonTrackLinks> >      t_iterL3OI_;
  edm::EDGetTokenT< std::vector<reco::MuonTrackLinks> >      t_iterL3IOFromL2_;
  edm::EDGetTokenT< std::vector<reco::MuonTrackLinks> >      t_iterL3FromL2_;
  edm::EDGetTokenT< std::vector<reco::Track> >               t_iterL3IOFromL1_;
  edm::EDGetTokenT< std::vector<reco::Muon> >                t_iterL3MuonNoID_;
  edm::EDGetTokenT< std::vector<reco::Muon> >                t_iterL3Muon_;

  edm::EDGetTokenT< reco::VertexCollection >                 t_hltIterL3MuonTrimmedPixelVertices_;
  edm::EDGetTokenT< reco::VertexCollection >                 t_hltIterL3FromL1MuonTrimmedPixelVertices_;

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

  edm::EDGetTokenT< LumiScalersCollection >                  t_lumiScaler_;
  edm::EDGetTokenT< LumiScalersCollection >                  t_offlineLumiScaler_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >         t_PUSummaryInfo_;
  edm::EDGetTokenT< GenEventInfoProduct >                    t_genEventInfo_;
  edm::EDGetTokenT< reco::GenParticleCollection >            t_genParticle_;

  TTree *ntuple_;
  static const int arrSize_ = 5000;

  // -- general event information
  bool isRealData_;
  int runNum_;
  int lumiBlockNum_;
  unsigned long long eventNum_;

  int nVertex_;

  double bunchID_;
  double instLumi_;
  double dataPU_;
  double dataPURMS_;
  double bunchLumi_;
  double offlineInstLumi_;
  double offlineDataPU_;
  double offlineDataPURMS_;
  double offlineBunchLumi_;
  int truePU_;
  double genEventWeight_;

  // -- generator level particles (only MC)
  int nGenParticle_;
  int genParticle_ID_[arrSize_];
  int genParticle_status_[arrSize_];
  int genParticle_mother_[arrSize_];

  double genParticle_pt_[arrSize_];
  double genParticle_eta_[arrSize_];
  double genParticle_phi_[arrSize_];
  double genParticle_px_[arrSize_];
  double genParticle_py_[arrSize_];
  double genParticle_pz_[arrSize_];
  double genParticle_energy_[arrSize_];
  double genParticle_charge_[arrSize_];

  int genParticle_isPrompt_[arrSize_];
  int genParticle_isPromptFinalState_[arrSize_];
  int genParticle_isTauDecayProduct_[arrSize_];
  int genParticle_isPromptTauDecayProduct_[arrSize_];
  int genParticle_isDirectPromptTauDecayProductFinalState_[arrSize_];
  int genParticle_isHardProcess_[arrSize_];
  int genParticle_isLastCopy_[arrSize_];
  int genParticle_isLastCopyBeforeFSR_[arrSize_];
  int genParticle_isPromptDecayed_[arrSize_];
  int genParticle_isDecayedLeptonHadron_[arrSize_];
  int genParticle_fromHardProcessBeforeFSR_[arrSize_];
  int genParticle_fromHardProcessDecayed_[arrSize_];
  int genParticle_fromHardProcessFinalState_[arrSize_];
  int genParticle_isMostlyLikePythia6Status3_[arrSize_];

  // -- trigger info.
  vector< std::string > vec_firedTrigger_;
  vector< std::string > vec_filterName_;
  vector< double > vec_HLTObj_pt_;
  vector< double > vec_HLTObj_eta_;
  vector< double > vec_HLTObj_phi_;

  vector< std::string > vec_myFiredTrigger_;
  vector< std::string > vec_myFilterName_;
  vector< double > vec_myHLTObj_pt_;
  vector< double > vec_myHLTObj_eta_;
  vector< double > vec_myHLTObj_phi_;

  // all L1 tracks
  vector<float>* m_trk_pt;
  vector<float>* m_trk_eta;
  vector<float>* m_trk_phi;
  vector<float>* m_trk_d0;   // (filled if L1Tk_nPar==5, else 999)
  vector<float>* m_trk_z0;
  vector<float>* m_trk_chi2;
  vector<float>* m_trk_bendchi2;
  vector<int>*   m_trk_nstub;
  vector<int>*   m_trk_lhits;
  vector<int>*   m_trk_dhits;
  vector<int>*   m_trk_seed;
  vector<unsigned int>*   m_trk_phiSector;
  vector<int>*   m_trk_genuine;
  vector<int>*   m_trk_loose;
  vector<int>*   m_trk_unknown;
  vector<int>*   m_trk_combinatoric;
  vector<int>*   m_trk_fake; //0 fake, 1 track from primary interaction, 2 secondary track
  vector<int>*   m_trk_matchtp_pdgid;
  vector<float>* m_trk_matchtp_pt;
  vector<float>* m_trk_matchtp_eta;
  vector<float>* m_trk_matchtp_phi;
  vector<float>* m_trk_matchtp_z0;
  vector<float>* m_trk_matchtp_dxy;

  vector<float>* m_stub_x;
  vector<float>* m_stub_y;
  vector<float>* m_stub_z;
  vector<int>*   m_stub_isBarrel; // stub is in barrel (1) or in disk (0)
  vector<int>*   m_stub_layer;

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

  // -- offline muon
  int nMuon_;

  double muon_pt_[arrSize_];
  double muon_eta_[arrSize_];
  double muon_phi_[arrSize_];
  double muon_px_[arrSize_];
  double muon_py_[arrSize_];
  double muon_pz_[arrSize_];
  double muon_dB_[arrSize_];
  double muon_charge_[arrSize_];
  int muon_isGLB_[arrSize_];
  int muon_isSTA_[arrSize_];
  int muon_isTRK_[arrSize_];
  int muon_isPF_[arrSize_];
  int muon_isTight_[arrSize_];
  int muon_isMedium_[arrSize_];
  int muon_isLoose_[arrSize_];
  int muon_isHighPt_[arrSize_];
  int muon_isHighPtNew_[arrSize_];
  int muon_isSoft_[arrSize_];

  double muon_iso03_sumPt_[arrSize_];
  double muon_iso03_hadEt_[arrSize_];
  double muon_iso03_emEt_[arrSize_];

  double muon_PFIso03_charged_[arrSize_];
  double muon_PFIso03_neutral_[arrSize_];
  double muon_PFIso03_photon_[arrSize_];
  double muon_PFIso03_sumPU_[arrSize_];

  double muon_PFIso04_charged_[arrSize_];
  double muon_PFIso04_neutral_[arrSize_];
  double muon_PFIso04_photon_[arrSize_];
  double muon_PFIso04_sumPU_[arrSize_];

  double muon_PFCluster03_ECAL_[arrSize_];
  double muon_PFCluster03_HCAL_[arrSize_];

  double muon_PFCluster04_ECAL_[arrSize_];
  double muon_PFCluster04_HCAL_[arrSize_];

  double muon_normChi2_global_[arrSize_];
  int muon_nTrackerHit_global_[arrSize_];
  int muon_nTrackerLayer_global_[arrSize_];
  int muon_nPixelHit_global_[arrSize_];
  int muon_nMuonHit_global_[arrSize_];

  double muon_normChi2_inner_[arrSize_];
  int muon_nTrackerHit_inner_[arrSize_];
  int muon_nTrackerLayer_inner_[arrSize_];
  int muon_nPixelHit_inner_[arrSize_];

  double muon_pt_tuneP_[arrSize_];
  double muon_ptError_tuneP_[arrSize_];

  double muon_dxyVTX_best_[arrSize_];
  double muon_dzVTX_best_[arrSize_];

  int muon_nMatchedStation_[arrSize_];
  int muon_nMatchedRPCLayer_[arrSize_];
  int muon_stationMask_[arrSize_];

  std::map<tmpTSOD,unsigned int> MuonIterSeedMap;
  std::map<tmpTSOD,unsigned int> hltIterL3OIMuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter0IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter2IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter3IterL3MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter0IterL3FromL1MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter2IterL3FromL1MuonTrackMap;
  std::map<tmpTSOD,unsigned int> hltIter3IterL3FromL1MuonTrackMap;

  // -- L3 muon
  int nL3Muon_;
  double L3Muon_pt_[arrSize_];
  double L3Muon_eta_[arrSize_];
  double L3Muon_phi_[arrSize_];
  double L3Muon_charge_[arrSize_];
  double L3Muon_trkPt_[arrSize_];

  // -- L2 muon
  int nL2Muon_;
  double L2Muon_pt_[arrSize_];
  double L2Muon_eta_[arrSize_];
  double L2Muon_phi_[arrSize_];
  double L2Muon_charge_[arrSize_];
  double L2Muon_trkPt_[arrSize_];

  // -- L1 muon
  int nL1Muon_;
  double L1Muon_pt_[arrSize_];
  double L1Muon_eta_[arrSize_];
  double L1Muon_phi_[arrSize_];
  double L1Muon_charge_[arrSize_];
  double L1Muon_quality_[arrSize_];

  // -- Tracker muon
  int nTkMuon_;
  double TkMuon_pt_[arrSize_];
  double TkMuon_eta_[arrSize_];
  double TkMuon_phi_[arrSize_];
  double TkMuon_charge_[arrSize_];
  double TkMuon_trkPt_[arrSize_];

  int    nIterL3OI_;
  double iterL3OI_inner_pt_[arrSize_];
  double iterL3OI_inner_eta_[arrSize_];
  double iterL3OI_inner_phi_[arrSize_];
  double iterL3OI_inner_charge_[arrSize_];
  double iterL3OI_outer_pt_[arrSize_];
  double iterL3OI_outer_eta_[arrSize_];
  double iterL3OI_outer_phi_[arrSize_];
  double iterL3OI_outer_charge_[arrSize_];
  double iterL3OI_global_pt_[arrSize_];
  double iterL3OI_global_eta_[arrSize_];
  double iterL3OI_global_phi_[arrSize_];
  double iterL3OI_global_charge_[arrSize_];

  int    nIterL3IOFromL2_;
  double iterL3IOFromL2_inner_pt_[arrSize_];
  double iterL3IOFromL2_inner_eta_[arrSize_];
  double iterL3IOFromL2_inner_phi_[arrSize_];
  double iterL3IOFromL2_inner_charge_[arrSize_];
  double iterL3IOFromL2_outer_pt_[arrSize_];
  double iterL3IOFromL2_outer_eta_[arrSize_];
  double iterL3IOFromL2_outer_phi_[arrSize_];
  double iterL3IOFromL2_outer_charge_[arrSize_];
  double iterL3IOFromL2_global_pt_[arrSize_];
  double iterL3IOFromL2_global_eta_[arrSize_];
  double iterL3IOFromL2_global_phi_[arrSize_];
  double iterL3IOFromL2_global_charge_[arrSize_];

  // -- iterL3 object from outside-in + inside-out step (from L2)
  int    nIterL3FromL2_;
  double iterL3FromL2_inner_pt_[arrSize_];
  double iterL3FromL2_inner_eta_[arrSize_];
  double iterL3FromL2_inner_phi_[arrSize_];
  double iterL3FromL2_inner_charge_[arrSize_];
  double iterL3FromL2_outer_pt_[arrSize_];
  double iterL3FromL2_outer_eta_[arrSize_];
  double iterL3FromL2_outer_phi_[arrSize_];
  double iterL3FromL2_outer_charge_[arrSize_];
  double iterL3FromL2_global_pt_[arrSize_];
  double iterL3FromL2_global_eta_[arrSize_];
  double iterL3FromL2_global_phi_[arrSize_];
  double iterL3FromL2_global_charge_[arrSize_];

  int    nIterL3IOFromL1_;
  double iterL3IOFromL1_pt_[arrSize_];
  double iterL3IOFromL1_eta_[arrSize_];
  double iterL3IOFromL1_phi_[arrSize_];
  double iterL3IOFromL1_charge_[arrSize_];

  // -- iterL3 object before applying ID @ HLT
  int nIterL3MuonNoID_;
  double iterL3MuonNoID_pt_[arrSize_];
  double iterL3MuonNoID_eta_[arrSize_];
  double iterL3MuonNoID_phi_[arrSize_];
  double iterL3MuonNoID_charge_[arrSize_];
  int iterL3MuonNoID_isGLB_[arrSize_];
  int iterL3MuonNoID_isSTA_[arrSize_];
  int iterL3MuonNoID_isTRK_[arrSize_];

  // -- iterL3 object after applying ID @ HLT
  int nIterL3Muon_;
  double iterL3Muon_pt_[arrSize_];
  double iterL3Muon_innerPt_[arrSize_];
  double iterL3Muon_eta_[arrSize_];
  double iterL3Muon_phi_[arrSize_];
  double iterL3Muon_charge_[arrSize_];
  int iterL3Muon_isGLB_[arrSize_];
  int iterL3Muon_isSTA_[arrSize_];
  int iterL3Muon_isTRK_[arrSize_];

  class seedTemplate {
  private:
    int nSeeds_;
    std::vector<int> dir_;
    std::vector<uint32_t> tsos_detId_;
    std::vector<float> tsos_pt_;
    std::vector<float> tsos_pt_val_;
    std::vector<float> tsos_eta_;
    std::vector<float> tsos_phi_;
    std::vector<float> tsos_glob_x_;
    std::vector<float> tsos_glob_y_;
    std::vector<float> tsos_glob_z_;
    std::vector<int> tsos_hasErr_;
    std::vector<float> tsos_err0_;
    std::vector<float> tsos_err1_;
    std::vector<float> tsos_err2_;
    std::vector<float> tsos_err3_;
    std::vector<float> tsos_err4_;
    std::vector<float> tsos_err5_;
    std::vector<float> tsos_err6_;
    std::vector<float> tsos_err7_;
    std::vector<float> tsos_err8_;
    std::vector<float> tsos_err9_;
    std::vector<float> tsos_err10_;
    std::vector<float> tsos_err11_;
    std::vector<float> tsos_err12_;
    std::vector<float> tsos_err13_;
    std::vector<float> tsos_err14_;
    std::vector<float> tsos_x_;
    std::vector<float> tsos_y_;
    std::vector<float> tsos_dxdz_;
    std::vector<float> tsos_dydz_;
    std::vector<float> tsos_px_;
    std::vector<float> tsos_py_;
    std::vector<float> tsos_pz_;
    std::vector<float> tsos_qbp_;
    std::vector<int> tsos_charge_;
    std::vector<int> iterL3Matched_;
    std::vector<int> iterL3Ref_;
    std::vector<int> tmpL3Ref_;
  public:
    void clear() {
      nSeeds_ = 0;
      dir_.clear();
      tsos_detId_.clear();
      tsos_pt_.clear();
      tsos_pt_val_.clear();
      tsos_eta_.clear();
      tsos_phi_.clear();
      tsos_glob_x_.clear();
      tsos_glob_y_.clear();
      tsos_glob_z_.clear();
      tsos_hasErr_.clear();
      tsos_err0_.clear();
      tsos_err1_.clear();
      tsos_err2_.clear();
      tsos_err3_.clear();
      tsos_err4_.clear();
      tsos_err5_.clear();
      tsos_err6_.clear();
      tsos_err7_.clear();
      tsos_err8_.clear();
      tsos_err9_.clear();
      tsos_err10_.clear();
      tsos_err11_.clear();
      tsos_err12_.clear();
      tsos_err13_.clear();
      tsos_err14_.clear();
      tsos_x_.clear();
      tsos_y_.clear();
      tsos_dxdz_.clear();
      tsos_dydz_.clear();
      tsos_px_.clear();
      tsos_py_.clear();
      tsos_pz_.clear();
      tsos_qbp_.clear();
      tsos_charge_.clear();
      iterL3Matched_.clear();
      iterL3Ref_.clear();
      tmpL3Ref_.clear();

      return;
    }

    void setBranch(TTree* tmpntpl, TString name) {
      tmpntpl->Branch("n"+name,             &nSeeds_);
      tmpntpl->Branch(name+"_dir",          &dir_);
      tmpntpl->Branch(name+"_tsos_detId",   &tsos_detId_);
      tmpntpl->Branch(name+"_tsos_pt",      &tsos_pt_);
      tmpntpl->Branch(name+"_tsos_pt_val",      &tsos_pt_val_);
      tmpntpl->Branch(name+"_tsos_eta",      &tsos_eta_);
      tmpntpl->Branch(name+"_tsos_phi",      &tsos_phi_);
      tmpntpl->Branch(name+"_tsos_glob_x",      &tsos_glob_x_);
      tmpntpl->Branch(name+"_tsos_glob_y",      &tsos_glob_y_);
      tmpntpl->Branch(name+"_tsos_glob_z",      &tsos_glob_z_);
      tmpntpl->Branch(name+"_tsos_hasErr",  &tsos_hasErr_);
      tmpntpl->Branch(name+"_tsos_err0",    &tsos_err0_);
      tmpntpl->Branch(name+"_tsos_err1",    &tsos_err1_);
      tmpntpl->Branch(name+"_tsos_err2",    &tsos_err2_);
      tmpntpl->Branch(name+"_tsos_err3",    &tsos_err3_);
      tmpntpl->Branch(name+"_tsos_err4",    &tsos_err4_);
      tmpntpl->Branch(name+"_tsos_err5",    &tsos_err5_);
      tmpntpl->Branch(name+"_tsos_err6",    &tsos_err6_);
      tmpntpl->Branch(name+"_tsos_err7",    &tsos_err7_);
      tmpntpl->Branch(name+"_tsos_err8",    &tsos_err8_);
      tmpntpl->Branch(name+"_tsos_err9",    &tsos_err9_);
      tmpntpl->Branch(name+"_tsos_err10",   &tsos_err10_);
      tmpntpl->Branch(name+"_tsos_err11",   &tsos_err11_);
      tmpntpl->Branch(name+"_tsos_err12",   &tsos_err12_);
      tmpntpl->Branch(name+"_tsos_err13",   &tsos_err13_);
      tmpntpl->Branch(name+"_tsos_err14",   &tsos_err14_);
      tmpntpl->Branch(name+"_tsos_x",       &tsos_x_);
      tmpntpl->Branch(name+"_tsos_y",       &tsos_y_);
      tmpntpl->Branch(name+"_tsos_dxdz",    &tsos_dxdz_);
      tmpntpl->Branch(name+"_tsos_dydz",    &tsos_dydz_);
      tmpntpl->Branch(name+"_tsos_px",      &tsos_px_);
      tmpntpl->Branch(name+"_tsos_py",      &tsos_py_);
      tmpntpl->Branch(name+"_tsos_pz",      &tsos_pz_);
      tmpntpl->Branch(name+"_tsos_qbp",     &tsos_qbp_);
      tmpntpl->Branch(name+"_tsos_charge",  &tsos_charge_);
      tmpntpl->Branch(name+"_iterL3Matched", &iterL3Matched_);
      tmpntpl->Branch(name+"_iterL3Ref", &iterL3Ref_);
      tmpntpl->Branch(name+"_tmpL3Ref", &tmpL3Ref_);

      return;
    }

    void fill(TrajectorySeed seed, edm::ESHandle<TrackerGeometry> tracker) {
      GlobalVector p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
      GlobalPoint x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

      dir_.push_back(seed.direction());
      tsos_detId_.push_back(seed.startingState().detId());
      tsos_pt_.push_back(seed.startingState().pt());
      tsos_pt_val_.push_back(p.perp());
      tsos_eta_.push_back(p.eta());
      tsos_phi_.push_back(p.phi());
      tsos_glob_x_.push_back(x.x());
      tsos_glob_y_.push_back(x.y());
      tsos_glob_z_.push_back(x.z());
      tsos_hasErr_.push_back(seed.startingState().hasError());
      if( seed.startingState().hasError() ) {
        tsos_err0_.push_back(seed.startingState().error(0));
        tsos_err1_.push_back(seed.startingState().error(1));
        tsos_err2_.push_back(seed.startingState().error(2));
        tsos_err3_.push_back(seed.startingState().error(3));
        tsos_err4_.push_back(seed.startingState().error(4));
        tsos_err5_.push_back(seed.startingState().error(5));
        tsos_err6_.push_back(seed.startingState().error(6));
        tsos_err7_.push_back(seed.startingState().error(7));
        tsos_err8_.push_back(seed.startingState().error(8));
        tsos_err9_.push_back(seed.startingState().error(9));
        tsos_err10_.push_back(seed.startingState().error(10));
        tsos_err11_.push_back(seed.startingState().error(11));
        tsos_err12_.push_back(seed.startingState().error(12));
        tsos_err13_.push_back(seed.startingState().error(13));
        tsos_err14_.push_back(seed.startingState().error(14));
      }
      else {  // sync vector size
        tsos_err0_.push_back(-99999.);
        tsos_err1_.push_back(-99999.);
        tsos_err2_.push_back(-99999.);
        tsos_err3_.push_back(-99999.);
        tsos_err4_.push_back(-99999.);
        tsos_err5_.push_back(-99999.);
        tsos_err6_.push_back(-99999.);
        tsos_err7_.push_back(-99999.);
        tsos_err8_.push_back(-99999.);
        tsos_err9_.push_back(-99999.);
        tsos_err10_.push_back(-99999.);
        tsos_err11_.push_back(-99999.);
        tsos_err12_.push_back(-99999.);
        tsos_err13_.push_back(-99999.);
        tsos_err14_.push_back(-99999.);
      }
      tsos_x_.push_back(seed.startingState().parameters().position().x());
      tsos_y_.push_back(seed.startingState().parameters().position().y());
      tsos_dxdz_.push_back(seed.startingState().parameters().dxdz());
      tsos_dydz_.push_back(seed.startingState().parameters().dydz());
      tsos_px_.push_back(seed.startingState().parameters().momentum().x());
      tsos_py_.push_back(seed.startingState().parameters().momentum().y());
      tsos_pz_.push_back(seed.startingState().parameters().momentum().z());
      tsos_qbp_.push_back(seed.startingState().parameters().qbp());
      tsos_charge_.push_back(seed.startingState().parameters().charge());
      nSeeds_++;

      return;
    }

    void linkIterL3(int linkNo) { iterL3Ref_.push_back(linkNo); }
    void linktmpL3(int linkNo) { tmpL3Ref_.push_back(linkNo); }
    void matchingL3(int linkNo) { iterL3Matched_.push_back(linkNo); }
  };

  class tmpTrk {
  private:
    double trkPt;
    double trkEta;
    double trkPhi;
    int trkCharge;
  public:
    bool isMatched(const reco::Track trk_) {
      return std::abs(trk_.pt()-trkPt)/trkPt < 0.01 && std::abs(trk_.eta()-trkEta) < 0.01 && std::abs(trk_.phi()-trkPhi) < 0.01 && trk_.charge()==trkCharge;
    }

    void fill(const reco::TrackRef trk_) {
      trkPt = trk_->pt();
      trkEta = trk_->eta();
      trkPhi = trk_->phi();
      trkCharge = trk_->charge();
    }

    tmpTrk(const reco::TrackRef trk_) { fill(trk_); }
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
  };

  class vtxTemplate {
  private:
    int nVtxs;
    std::vector<int> isValid;
    std::vector<double> chi2;
    std::vector<double> ndof;
    std::vector<double> nTracks;
    std::vector<double> vtxX;
    std::vector<double> vtxXerr;
    std::vector<double> vtxY;
    std::vector<double> vtxYerr;
    std::vector<double> vtxZ;
    std::vector<double> vtxZerr;
  public:
    void clear() {
      nVtxs = 0;
      isValid.clear();
      chi2.clear();
      ndof.clear();
      nTracks.clear();
      vtxX.clear();
      vtxXerr.clear();
      vtxY.clear();
      vtxYerr.clear();
      vtxZ.clear();
      vtxZerr.clear();

      return;
    }

    void setBranch(TTree* tmpntpl, TString name) {
      tmpntpl->Branch("n"+name, &nVtxs);
      tmpntpl->Branch(name+"_isValid", &isValid);
      tmpntpl->Branch(name+"_chi2", &chi2);
      tmpntpl->Branch(name+"_ndof", &ndof);
      tmpntpl->Branch(name+"_nTracks", &nTracks);
      tmpntpl->Branch(name+"_x", &vtxX);
      tmpntpl->Branch(name+"_xerr", &vtxXerr);
      tmpntpl->Branch(name+"_y", &vtxY);
      tmpntpl->Branch(name+"_yerr", &vtxYerr);
      tmpntpl->Branch(name+"_z", &vtxZ);
      tmpntpl->Branch(name+"_zerr", &vtxZerr);

      return;
    }

    void fill(const reco::Vertex vtx) {
      isValid.push_back(vtx.isValid());
      chi2.push_back(vtx.chi2());
      ndof.push_back(vtx.ndof());
      nTracks.push_back(vtx.nTracks());
      vtxX.push_back(vtx.x());
      vtxXerr.push_back(vtx.xError());
      vtxY.push_back(vtx.y());
      vtxYerr.push_back(vtx.yError());
      vtxZ.push_back(vtx.z());
      vtxZerr.push_back(vtx.zError());
      nVtxs++;

      return;
    }
  };

  class hitTemplate {
  private:
    int nHits;
    std::vector<bool> isValid;
    std::vector<float> localx, localy, localz, globalx, globaly, globalz;
    // std::vector<float> localxx, localyy, localxy;
    // std::vector<float> globalcxx, globalcyx, globalcyy, globalczx, globalczy, globalczz, globalctx, globalcty, blobalctz, globalctt;
  public:
    void clear() {
      nHits = 0;
      isValid.clear();
      localx.clear(); localy.clear(); localz.clear(); globalx.clear(); globaly.clear(); globalz.clear();
      // localxx.clear(); localyy.clear(); localxy.clear();
      // globalcxx.clear(); globalcyx.clear(); globalcyy.clear(); globalczx.clear(); globalczy.clear();
      // globalczz.clear(); globalctx.clear(); globalcty.clear(); blobalctz.clear(); globalctt.clear();

      return;
    }

    void fill(int idx, TrackingRecHit* rechit) {
      isValid.push_back(rechit->isValid());
      localx.push_back(rechit->localPosition().x());
      localy.push_back(rechit->localPosition().y());
      localz.push_back(rechit->localPosition().z());
      globalx.push_back(rechit->globalPosition().x());
      globaly.push_back(rechit->globalPosition().y());
      globalz.push_back(rechit->globalPosition().z());

      return;
    }

    void setBranch(TTree* tmpntpl, TString name) {
      tmpntpl->Branch(name+"_isValid", &isValid);
      tmpntpl->Branch(name+"_localX", &localx);
      tmpntpl->Branch(name+"_localY", &localy);
      tmpntpl->Branch(name+"_localZ", &localz);
      tmpntpl->Branch(name+"_globalX", &globalx);
      tmpntpl->Branch(name+"_globalY", &globaly);
      tmpntpl->Branch(name+"_globalZ", &globalz);

      return;
    }
  };

  std::vector<tmpTrk> iterL3IDpassed;

  seedTemplate* SThltIterL3OISeedsFromL2Muons = new seedTemplate();
  seedTemplate* SThltIter0IterL3MuonPixelSeedsFromPixelTracks = new seedTemplate();
  seedTemplate* SThltIter2IterL3MuonPixelSeeds = new seedTemplate();
  seedTemplate* SThltIter3IterL3MuonPixelSeeds = new seedTemplate();
  seedTemplate* SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = new seedTemplate();
  seedTemplate* SThltIter2IterL3FromL1MuonPixelSeeds = new seedTemplate();
  seedTemplate* SThltIter3IterL3FromL1MuonPixelSeeds = new seedTemplate();

  trkTemplate* TThltIterL3OIMuonTrack = new trkTemplate();
  trkTemplate* TThltIter0IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter2IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter3IterL3MuonTrack = new trkTemplate();
  trkTemplate* TThltIter0IterL3FromL1MuonTrack = new trkTemplate();
  trkTemplate* TThltIter2IterL3FromL1MuonTrack = new trkTemplate();
  trkTemplate* TThltIter3IterL3FromL1MuonTrack = new trkTemplate();

  vtxTemplate* VThltIterL3MuonTrimmedPixelVertices = new vtxTemplate();
  vtxTemplate* VThltIterL3FromL1MuonTrimmedPixelVertices = new vtxTemplate();

  void fill_trackTemplate(const edm::Event &iEvent, edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
    edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_, edm::Handle<TrackingParticleCollection>& TPCollection_,
    std::map<tmpTSOD,unsigned int>& trkMap, trkTemplate* TTtrack);
};
