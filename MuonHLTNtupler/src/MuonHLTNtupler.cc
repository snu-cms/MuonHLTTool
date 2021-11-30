// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTNtupler.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/Utils/interface/associationMapFilterValues.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"




#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

using namespace std;
using namespace reco;
using namespace edm;


MuonHLTNtupler::MuonHLTNtupler(const edm::ParameterSet& iConfig):
doMVA(iConfig.getParameter<bool>("doMVA")),
doSeed(iConfig.getParameter<bool>("doSeed")),
DebugMode(iConfig.getParameter<bool>("DebugMode")),

associatorToken(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("associator"))),
trackingParticleToken(consumes<TrackingParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingParticle"))),
t_beamSpot_          ( consumes< reco::BeamSpot >                         (iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"     )) ),
t_offlineMuon_       ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
t_offlineVertex_     ( consumes< reco::VertexCollection >                 (iConfig.getUntrackedParameter<edm::InputTag>("offlineVertex"     )) ),
t_triggerResults_    ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults"    )) ),
t_triggerEvent_      ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"      )) ),
t_myTriggerResults_  ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerResults"  )) ),
t_myTriggerEvent_    ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerEvent"    )) ),
t_L3Muon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L3Muon"            )) ),
t_L2Muon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2Muon"            )) ),
t_L1Muon_            ( consumes< l1t::MuonBxCollection  >                 (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon"            )) ),
t_TkMuon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("TkMuon"            )) ),

t_iterL3OI_          ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3OI"          )) ),
t_iterL3IOFromL2_    ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3IOFromL2"    )) ),
t_iterL3FromL2_      ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3FromL2"      )) ),
t_iterL3IOFromL1_    ( consumes< std::vector<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("iterL3IOFromL1"    )) ),
t_iterL3MuonNoID_    ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("iterL3MuonNoID"    )) ),
t_iterL3Muon_        ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("iterL3Muon"        )) ),

t_hltIterL3MuonTrimmedPixelVertices_       ( consumes< reco::VertexCollection > (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3MuonTrimmedPixelVertices"     )) ),
t_hltIterL3FromL1MuonTrimmedPixelVertices_ ( consumes< reco::VertexCollection > (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3FromL1MuonTrimmedPixelVertices"     )) ),

t_hltIterL3OISeedsFromL2Muons_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3OISeedsFromL2Muons")) ),
t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3MuonPixelSeedsFromPixelTracks")) ),
t_hltIter2IterL3MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3MuonPixelSeeds")) ),
t_hltIter3IterL3MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3MuonPixelSeeds")) ),
t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks")) ),
t_hltIter2IterL3FromL1MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3FromL1MuonPixelSeeds")) ),
t_hltIter3IterL3FromL1MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3FromL1MuonPixelSeeds")) ),

t_hltIterL3OIMuonTrack_    ( consumes< edm::View<reco::Track> >                  (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3OIMuonTrack"    )) ),
t_hltIter0IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3MuonTrack"    )) ),
t_hltIter2IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3MuonTrack"    )) ),
t_hltIter3IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3MuonTrack"    )) ),
t_hltIter0IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3FromL1MuonTrack"    )) ),
t_hltIter2IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3FromL1MuonTrack"    )) ),
t_hltIter3IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3FromL1MuonTrack"    )) ),

t_lumiScaler_        ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("lumiScaler"        )) ),
t_offlineLumiScaler_ ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("offlineLumiScaler" )) ),
t_PUSummaryInfo_     ( consumes< std::vector<PileupSummaryInfo> >         (iConfig.getUntrackedParameter<edm::InputTag>("PUSummaryInfo"     )) ),
t_genEventInfo_      ( consumes< GenEventInfoProduct >                    (iConfig.getUntrackedParameter<edm::InputTag>("genEventInfo"      )) ),
t_genParticle_       ( consumes< edm::View<reco::GenParticle> >           (iConfig.getUntrackedParameter<edm::InputTag>("genParticle"       )) ),

t_l1Matches_           ( consumes< pat::TriggerObjectStandAloneMatch >                (iConfig.getParameter<edm::InputTag>("l1Matches"))),
t_l1MatchesQuality_    ( consumes< edm::ValueMap<int> >                               (iConfig.getParameter<edm::InputTag>("l1MatchesQuality"))),
t_l1MatchesDeltaR_     ( consumes< edm::ValueMap<float> >                             (iConfig.getParameter<edm::InputTag>("l1MatchesDeltaR"))),
t_l1MatchesByQ_        ( consumes< pat::TriggerObjectStandAloneMatch >                (iConfig.getParameter<edm::InputTag>("l1MatchesByQ"))),
t_l1MatchesByQQuality_ ( consumes< edm::ValueMap<int> >                               (iConfig.getParameter<edm::InputTag>("l1MatchesByQQuality"))),
t_l1MatchesByQDeltaR_  ( consumes< edm::ValueMap<float> >                             (iConfig.getParameter<edm::InputTag>("l1MatchesByQDeltaR"))),
bs(0)
{
  trackCollectionNames_   = iConfig.getUntrackedParameter<std::vector<std::string>   >("trackCollectionNames");
  trackCollectionLabels_  = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("trackCollectionLabels");
  associationLabels_      = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("associationLabels");
  if( trackCollectionNames_.size() != trackCollectionLabels_.size() || trackCollectionLabels_.size() != associationLabels_.size()) {
    throw cms::Exception("ConfigurationError")
      << "Number of track collection names is different from number of track collection names or association labels";
  }
  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    trackCollectionTokens_.push_back(     consumes< edm::View<reco::Track> >(  trackCollectionLabels_[i]) );
    simToRecoCollectionTokens_.push_back( consumes<reco::SimToRecoCollection>( associationLabels_[i]) );
    recoToSimCollectionTokens_.push_back( consumes<reco::RecoToSimCollection>( associationLabels_[i]) );
    trkTemplates_.push_back( new trkTemplate() );
    tpTemplates_.push_back(  new tpTemplate()  );
  }

  mvaFileHltIter2IterL3MuonPixelSeeds_B_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_B");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaFileHltIter2IterL3MuonPixelSeeds_E_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_E");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E");

  mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3MuonPixelSeeds_B");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E");
  mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3MuonPixelSeeds_E");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E");
}

void MuonHLTNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init();

  // -- basic info.
  isRealData_ = iEvent.isRealData();

  runNum_       = iEvent.id().run();
  lumiBlockNum_ = iEvent.id().luminosityBlock();
  eventNum_     = iEvent.id().event();

  edm::Handle<reco::BeamSpot> h_beamSpot;
  iEvent.getByToken(t_beamSpot_, h_beamSpot);
  bs = h_beamSpot.isValid() ? &*h_beamSpot : 0;
  bs_x0_ = bs->x0();
  bs_y0_ = bs->y0();
  bs_z0_ = bs->z0();
  bs_sigmaZ_ = bs->sigmaZ();
  bs_dxdz_ = bs->dxdz();
  bs_dydz_ = bs->dydz();
  bs_x0Error_ = bs->x0Error();
  bs_y0Error_ = bs->y0Error();
  bs_z0Error_ = bs->z0Error();
  bs_sigmaZ0Error_ = bs->sigmaZ0Error();
  bs_dxdzError_ = bs->dxdzError();
  bs_dydzError_ = bs->dydzError();


  // -- vertex
  edm::Handle<reco::VertexCollection> h_offlineVertex;
  if( iEvent.getByToken(t_offlineVertex_, h_offlineVertex) )
  {
    int nGoodVtx = 0;
    for(reco::VertexCollection::const_iterator it = h_offlineVertex->begin(); it != h_offlineVertex->end(); ++it)
      if( it->isValid() ) nGoodVtx++;

    nVertex_ = nGoodVtx;
  }

  // -- hltIterL3MuonTrimmedPixelVertices
  edm::Handle<reco::VertexCollection> h_hltIterL3MuonTrimmedPixelVertices;
  if( iEvent.getByToken(t_hltIterL3MuonTrimmedPixelVertices_, h_hltIterL3MuonTrimmedPixelVertices) )
  {
    for(reco::VertexCollection::const_iterator it = h_hltIterL3MuonTrimmedPixelVertices->begin(); it != h_hltIterL3MuonTrimmedPixelVertices->end(); ++it)
      VThltIterL3MuonTrimmedPixelVertices->fill(*it);
  }

  // -- hltIterL3FromL1MuonTrimmedPixelVertices
  edm::Handle<reco::VertexCollection> h_hltIterL3FromL1MuonTrimmedPixelVertices;
  if( iEvent.getByToken(t_hltIterL3FromL1MuonTrimmedPixelVertices_, h_hltIterL3FromL1MuonTrimmedPixelVertices) )
  {
    for(reco::VertexCollection::const_iterator it = h_hltIterL3FromL1MuonTrimmedPixelVertices->begin(); it != h_hltIterL3FromL1MuonTrimmedPixelVertices->end(); ++it)
      VThltIterL3FromL1MuonTrimmedPixelVertices->fill(*it);
  }

  if( isRealData_ )
  {
    bunchID_ = iEvent.bunchCrossing();

    // -- lumi scaler @ HLT
    edm::Handle<LumiScalersCollection> h_lumiScaler;
    if( iEvent.getByToken(t_lumiScaler_, h_lumiScaler) && h_lumiScaler->begin() != h_lumiScaler->end() )
    {
      instLumi_  = h_lumiScaler->begin()->instantLumi();
      dataPU_    = h_lumiScaler->begin()->pileup();
      dataPURMS_ = h_lumiScaler->begin()->pileupRMS();
      bunchLumi_ = h_lumiScaler->begin()->bunchLumi();
    }

    // -- lumi scaler @ offline
    edm::Handle<LumiScalersCollection> h_offlineLumiScaler;
    if( iEvent.getByToken(t_offlineLumiScaler_, h_offlineLumiScaler) && h_offlineLumiScaler->begin() != h_offlineLumiScaler->end() )
    {
      offlineInstLumi_  = h_offlineLumiScaler->begin()->instantLumi();
      offlineDataPU_    = h_offlineLumiScaler->begin()->pileup();
      offlineDataPURMS_ = h_offlineLumiScaler->begin()->pileupRMS();
      offlineBunchLumi_ = h_offlineLumiScaler->begin()->bunchLumi();
    }
  }

  // -- True PU info: only for MC -- //
  if( !isRealData_ )
  {
    edm::Handle<std::vector< PileupSummaryInfo > > h_PUSummaryInfo;

    if( iEvent.getByToken(t_PUSummaryInfo_,h_PUSummaryInfo) )
    {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = h_PUSummaryInfo->begin(); PVI != h_PUSummaryInfo->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0)
        {
          truePU_ = PVI->getTrueNumInteractions();
          continue;
        }
      } // -- end of PU iteration -- //
    } // -- end of if ( token exists )
  } // -- end of isMC -- //

  // -- fill each object
  // Fill_L1Track(iEvent, iSetup);
  Fill_Muon(iEvent);
  Fill_HLT(iEvent, 0); // -- original HLT objects saved in data taking
  Fill_HLT(iEvent, 1); // -- rerun objects
  Fill_HLTMuon(iEvent);
  Fill_L1Muon(iEvent);
  Fill_IterL3(iEvent, iSetup);
  //if( doSeed )  Fill_Seed(iEvent, iSetup);
  if( !isRealData_ ) {
    Fill_GenParticle(iEvent);
    Fill_TP(iEvent, TrkParticle);
  }

  ntuple_->Fill();
}

void MuonHLTNtupler::beginJob()
{
  mvaHltIter2IterL3MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_B_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_, false, 7),
               new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_E_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_, false, 7) )
  };
  mvaHltIter2IterL3FromL1MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_, true, 7),
               new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_, true, 7) )
  };

  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");

  Make_Branch();
}

void MuonHLTNtupler::Init()
{
  isRealData_ = false;

  runNum_       = -999;
  lumiBlockNum_ = -999;
  eventNum_     = 0;

  bunchID_ = -999;

  bs_x0_ = -999;
  bs_y0_ = -999;
  bs_z0_ = -999;
  bs_sigmaZ_ = -999;
  bs_dxdz_ = -999;
  bs_dydz_ = -999;
  bs_x0Error_ = -999;
  bs_y0Error_ = -999;
  bs_z0Error_ = -999;
  bs_sigmaZ0Error_ = -999;
  bs_dxdzError_ = -999;
  bs_dydzError_ = -999;

  nVertex_ = -999;

  instLumi_  = -999;
  dataPU_    = -999;
  dataPURMS_ = -999;
  bunchLumi_ = -999;

  offlineInstLumi_  = -999;
  offlineDataPU_    = -999;
  offlineDataPURMS_ = -999;
  offlineBunchLumi_ = -999;

  truePU_ = -999;

  genEventWeight_ = -999;

  nGenParticle_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    genParticle_ID_[i] = -999;
    genParticle_status_[i] = -999;
    genParticle_mother_[i] = -999;

    genParticle_pt_[i]     = -999;
    genParticle_eta_[i]    = -999;
    genParticle_phi_[i]    = -999;
    genParticle_px_[i]     = -999;
    genParticle_py_[i]     = -999;
    genParticle_pz_[i]     = -999;
    genParticle_energy_[i] = -999;
    genParticle_charge_[i] = -999;

    genParticle_isPrompt_[i] = 0;
    genParticle_isPromptFinalState_[i] = 0;
    genParticle_isTauDecayProduct_[i] = 0;
    genParticle_isPromptTauDecayProduct_[i] = 0;
    genParticle_isDirectPromptTauDecayProductFinalState_[i] = 0;
    genParticle_isHardProcess_[i] = 0;
    genParticle_isLastCopy_[i] = 0;
    genParticle_isLastCopyBeforeFSR_[i] = 0;
    genParticle_isPromptDecayed_[i] = 0;
    genParticle_isDecayedLeptonHadron_[i] = 0;
    genParticle_fromHardProcessBeforeFSR_[i] = 0;
    genParticle_fromHardProcessDecayed_[i] = 0;
    genParticle_fromHardProcessFinalState_[i] = 0;
    genParticle_isMostlyLikePythia6Status3_[i] = 0;

    // -- L1 matched gen particles -- //
    genParticle_l1pt_[i]        = -999;
    genParticle_l1eta_[i]       = -999;
    genParticle_l1phi_[i]       = -999;
    genParticle_l1charge_[i]    = -999;
    genParticle_l1q_[i]         = -999;
    genParticle_l1dr_[i]        = -999;
    genParticle_l1ptByQ_[i]     = -999;
    genParticle_l1etaByQ_[i]    = -999;
    genParticle_l1phiByQ_[i]    = -999;
    genParticle_l1chargeByQ_[i] = -999;
    genParticle_l1qByQ_[i]      = -999;
    genParticle_l1drByQ_[i]     = -999;
  }

  // -- original trigger objects -- //
  vec_firedTrigger_.clear();
  vec_filterName_.clear();
  vec_HLTObj_pt_.clear();
  vec_HLTObj_eta_.clear();
  vec_HLTObj_phi_.clear();

  // -- HLT rerun objects -- //
  vec_myFiredTrigger_.clear();
  vec_myFilterName_.clear();
  vec_myHLTObj_pt_.clear();
  vec_myHLTObj_eta_.clear();
  vec_myHLTObj_phi_.clear();

  MuonIterSeedMap.clear();
  hltIterL3OIMuonTrackMap.clear();
  hltIter0IterL3MuonTrackMap.clear();
  hltIter2IterL3MuonTrackMap.clear();
  hltIter3IterL3MuonTrackMap.clear();
  hltIter0IterL3FromL1MuonTrackMap.clear();
  hltIter2IterL3FromL1MuonTrackMap.clear();
  hltIter3IterL3FromL1MuonTrackMap.clear();
  iterL3IDpassed.clear();
  iterL3NoIDpassed.clear();

  nMuon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    muon_pt_[i] = -999;
    muon_eta_[i] = -999;
    muon_phi_[i] = -999;
    muon_px_[i] = -999;
    muon_py_[i] = -999;
    muon_pz_[i] = -999;
    muon_dB_[i] = -999;
    muon_charge_[i] = -999;
    muon_isGLB_[i] = 0;
    muon_isSTA_[i] = 0;
    muon_isTRK_[i] = 0;
    muon_isPF_[i] = 0;
    muon_isTight_[i] = 0;
    muon_isMedium_[i] = 0;
    muon_isLoose_[i] = 0;
    muon_isHighPt_[i] = 0;
    muon_isHighPtNew_[i] = 0;
    muon_isSoft_[i] = 0;

    muon_iso03_sumPt_[i] = -999;
    muon_iso03_hadEt_[i] = -999;
    muon_iso03_emEt_[i] = -999;

    muon_PFIso03_charged_[i] = -999;
    muon_PFIso03_neutral_[i] = -999;
    muon_PFIso03_photon_[i] = -999;
    muon_PFIso03_sumPU_[i] = -999;

    muon_PFIso04_charged_[i] = -999;
    muon_PFIso04_neutral_[i] = -999;
    muon_PFIso04_photon_[i] = -999;
    muon_PFIso04_sumPU_[i] = -999;

    muon_PFCluster03_ECAL_[i] = -999;
    muon_PFCluster03_HCAL_[i] = -999;

    muon_PFCluster04_ECAL_[i] = -999;
    muon_PFCluster04_HCAL_[i] = -999;

    muon_normChi2_global_[i] = -999;
    muon_nTrackerHit_global_[i] = -999;
    muon_nTrackerLayer_global_[i] = -999;
    muon_nPixelHit_global_[i] = -999;
    muon_nMuonHit_global_[i] = -999;

    muon_normChi2_inner_[i] = -999;
    muon_nTrackerHit_inner_[i] = -999;
    muon_nTrackerLayer_inner_[i] = -999;
    muon_nPixelHit_inner_[i] = -999;

    muon_pt_tuneP_[i] = -999;
    muon_ptError_tuneP_[i] = -999;

    muon_dxyVTX_best_[i] = -999;
    muon_dzVTX_best_[i] = -999;

    muon_nMatchedStation_[i] = -999;
    muon_nMatchedRPCLayer_[i] = -999;
    muon_stationMask_[i] = -999;

    muon_dxy_bs_[i] = -999;
    muon_dxyError_bs_[i] = -999;
    muon_dz_bs_[i] = -999;
    muon_dzError_[i] = -999;
    muon_IPSig_[i] = -999;
  }

  nL3Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L3Muon_pt_[i] = -999;
    L3Muon_eta_[i] = -999;
    L3Muon_phi_[i] = -999;
    L3Muon_charge_[i] = -999;
    L3Muon_trkPt_[i] = -999;
  }

  nL2Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L2Muon_pt_[i] = -999;
    L2Muon_eta_[i] = -999;
    L2Muon_phi_[i] = -999;
    L2Muon_charge_[i] = -999;
    L2Muon_trkPt_[i] = -999;
  }

  nTkMuon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    TkMuon_pt_[i] = -999;
    TkMuon_eta_[i] = -999;
    TkMuon_phi_[i] = -999;
    TkMuon_charge_[i] = -999;
    TkMuon_trkPt_[i] = -999;
  }

  nL1Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L1Muon_pt_[i] = -999;
    L1Muon_eta_[i] = -999;
    L1Muon_phi_[i] = -999;
    L1Muon_charge_[i] = -999;
    L1Muon_quality_[i] = -999;
    L1Muon_etaAtVtx_[i] = -999;
    L1Muon_phiAtVtx_[i] = -999;
  }

  nIterL3OI_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3OI_inner_pt_[i] = -999;
    iterL3OI_inner_eta_[i] = -999;
    iterL3OI_inner_phi_[i] = -999;
    iterL3OI_inner_charge_[i] = -999;
    iterL3OI_outer_pt_[i] = -999;
    iterL3OI_outer_eta_[i] = -999;
    iterL3OI_outer_phi_[i] = -999;
    iterL3OI_outer_charge_[i] = -999;
    iterL3OI_global_pt_[i] = -999;
    iterL3OI_global_eta_[i] = -999;
    iterL3OI_global_phi_[i] = -999;
    iterL3OI_global_charge_[i] = -999;
  }

  nIterL3IOFromL2_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3IOFromL2_inner_pt_[i] = -999;
    iterL3IOFromL2_inner_eta_[i] = -999;
    iterL3IOFromL2_inner_phi_[i] = -999;
    iterL3IOFromL2_inner_charge_[i] = -999;
    iterL3IOFromL2_outer_pt_[i] = -999;
    iterL3IOFromL2_outer_eta_[i] = -999;
    iterL3IOFromL2_outer_phi_[i] = -999;
    iterL3IOFromL2_outer_charge_[i] = -999;
    iterL3IOFromL2_global_pt_[i] = -999;
    iterL3IOFromL2_global_eta_[i] = -999;
    iterL3IOFromL2_global_phi_[i] = -999;
    iterL3IOFromL2_global_charge_[i] = -999;
  }

  nIterL3FromL2_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3FromL2_inner_pt_[i] = -999;
    iterL3FromL2_inner_eta_[i] = -999;
    iterL3FromL2_inner_phi_[i] = -999;
    iterL3FromL2_inner_charge_[i] = -999;
    iterL3FromL2_outer_pt_[i] = -999;
    iterL3FromL2_outer_eta_[i] = -999;
    iterL3FromL2_outer_phi_[i] = -999;
    iterL3FromL2_outer_charge_[i] = -999;
    iterL3FromL2_global_pt_[i] = -999;
    iterL3FromL2_global_eta_[i] = -999;
    iterL3FromL2_global_phi_[i] = -999;
    iterL3FromL2_global_charge_[i] = -999;
  }

  nIterL3IOFromL1_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3IOFromL1_pt_[i] = -999;
    iterL3IOFromL1_eta_[i] = -999;
    iterL3IOFromL1_phi_[i] = -999;
    iterL3IOFromL1_charge_[i] = -999;
  }


  nIterL3MuonNoID_ = 0;
  for (int i=0; i<arrSize_; ++i)
  {
    iterL3MuonNoID_pt_[i] = -999;
    iterL3MuonNoID_innerPt_[i] = -999;
    iterL3MuonNoID_eta_[i] = -999;
    iterL3MuonNoID_phi_[i] = -999;
    iterL3MuonNoID_charge_[i] = -999;

    iterL3MuonNoID_isGLB_[i] = 0;
    iterL3MuonNoID_isSTA_[i] = 0;
    iterL3MuonNoID_isTRK_[i] = 0;
  }

  nIterL3Muon_ = 0;
  for (int i=0; i<arrSize_; ++i)
  {
    iterL3Muon_pt_[i] = -999;
    iterL3Muon_innerPt_[i] = -999;
    iterL3Muon_eta_[i] = -999;
    iterL3Muon_phi_[i] = -999;
    iterL3Muon_charge_[i] = -999;

    iterL3Muon_isGLB_[i] = 0;
    iterL3Muon_isSTA_[i] = 0;
    iterL3Muon_isTRK_[i] = 0;
  }

  //SThltIterL3OISeedsFromL2Muons->clear();
  //SThltIter0IterL3MuonPixelSeedsFromPixelTracks->clear();
  //SThltIter2IterL3MuonPixelSeeds->clear();
  //SThltIter3IterL3MuonPixelSeeds->clear();
  //SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->clear();
  //SThltIter2IterL3FromL1MuonPixelSeeds->clear();
  //SThltIter3IterL3FromL1MuonPixelSeeds->clear();

  TThltIterL3OIMuonTrack->clear();
  TThltIter0IterL3MuonTrack->clear();
  TThltIter2IterL3MuonTrack->clear();
  TThltIter3IterL3MuonTrack->clear();
  TThltIter0IterL3FromL1MuonTrack->clear();
  TThltIter2IterL3FromL1MuonTrack->clear();
  TThltIter3IterL3FromL1MuonTrack->clear();

  TrkParticle->clear();

  VThltIterL3MuonTrimmedPixelVertices->clear();
  VThltIterL3FromL1MuonTrimmedPixelVertices->clear();

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    trkTemplates_.at(i)->clear();
    tpTemplates_.at(i)->clear();
  }

}

void MuonHLTNtupler::Make_Branch()
{
  ntuple_->Branch("isRealData", &isRealData_, "isRealData/O"); // -- O: boolean -- //
  ntuple_->Branch("runNum",&runNum_,"runNum/I");
  ntuple_->Branch("lumiBlockNum",&lumiBlockNum_,"lumiBlockNum/I");
  ntuple_->Branch("eventNum",&eventNum_,"eventNum/l"); // -- unsigned long long -- //
  ntuple_->Branch("bs_x0", &bs_x0_, "bs_x0/D");
  ntuple_->Branch("bs_y0", &bs_y0_, "bs_y0/D");
  ntuple_->Branch("bs_z0", &bs_z0_, "bs_z0/D");
  ntuple_->Branch("bs_sigmaZ", &bs_sigmaZ_, "bs_sigmaZ/D");
  ntuple_->Branch("bs_dxdz", &bs_dxdz_, "bs_dxdz/D");
  ntuple_->Branch("bs_dydz", &bs_dydz_, "bs_dydz/D");
  ntuple_->Branch("bs_x0Error", &bs_x0Error_, "bs_x0Error/D");
  ntuple_->Branch("bs_y0Error", &bs_y0Error_, "bs_y0Error/D");
  ntuple_->Branch("bs_z0Error", &bs_z0Error_, "bs_z0Error/D");
  ntuple_->Branch("bs_sigmaZ0Error", &bs_sigmaZ0Error_, "bs_sigmaZ0Error/D");
  ntuple_->Branch("bs_dxdzError", &bs_dxdzError_, "bs_dxdzError/D");
  ntuple_->Branch("bs_dydzError", &bs_dydzError_, "bs_dydzError/D");
  ntuple_->Branch("nVertex", &nVertex_, "nVertex/I");
  ntuple_->Branch("bunchID", &bunchID_, "bunchID/D");
  ntuple_->Branch("instLumi", &instLumi_, "instLumi/D");
  ntuple_->Branch("dataPU", &dataPU_, "dataPU/D");
  ntuple_->Branch("dataPURMS", &dataPURMS_, "dataPURMS/D");
  ntuple_->Branch("bunchLumi", &bunchLumi_, "bunchLumi/D");
  ntuple_->Branch("offlineInstLumi", &offlineInstLumi_, "offlineInstLumi/D");
  ntuple_->Branch("offlineDataPU", &offlineDataPU_, "offlineDataPU/D");
  ntuple_->Branch("offlineDataPURMS", &offlineDataPURMS_, "offlineDataPURMS/D");
  ntuple_->Branch("offlineBunchLumi", &offlineBunchLumi_, "offlineBunchLumi/D");
  ntuple_->Branch("truePU", &truePU_, "truePU/I");

  ntuple_->Branch("genEventWeight", &genEventWeight_, "genEventWeight/D");
  ntuple_->Branch("nGenParticle", &nGenParticle_, "nGenParticle/I");
  ntuple_->Branch("genParticle_ID", &genParticle_ID_, "genParticle_ID[nGenParticle]/I");
  ntuple_->Branch("genParticle_status", &genParticle_status_, "genParticle_status[nGenParticle]/I");
  ntuple_->Branch("genParticle_mother", &genParticle_mother_, "genParticle_mother[nGenParticle]/I");
  ntuple_->Branch("genParticle_pt", &genParticle_pt_, "genParticle_pt[nGenParticle]/D");
  ntuple_->Branch("genParticle_eta", &genParticle_eta_, "genParticle_eta[nGenParticle]/D");
  ntuple_->Branch("genParticle_phi", &genParticle_phi_, "genParticle_phi[nGenParticle]/D");
  ntuple_->Branch("genParticle_px", &genParticle_px_, "genParticle_px[nGenParticle]/D");
  ntuple_->Branch("genParticle_py", &genParticle_py_, "genParticle_py[nGenParticle]/D");
  ntuple_->Branch("genParticle_pz", &genParticle_pz_, "genParticle_pz[nGenParticle]/D");
  ntuple_->Branch("genParticle_energy", &genParticle_energy_, "genParticle_energy[nGenParticle]/D");
  ntuple_->Branch("genParticle_charge", &genParticle_charge_, "genParticle_charge[nGenParticle]/D");
  ntuple_->Branch("genParticle_isPrompt", &genParticle_isPrompt_, "genParticle_isPrompt[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptFinalState", &genParticle_isPromptFinalState_, "genParticle_isPromptFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isTauDecayProduct", &genParticle_isTauDecayProduct_, "genParticle_isTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptTauDecayProduct", &genParticle_isPromptTauDecayProduct_, "genParticle_isPromptTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDirectPromptTauDecayProductFinalState", &genParticle_isDirectPromptTauDecayProductFinalState_, "genParticle_isDirectPromptTauDecayProductFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isHardProcess", &genParticle_isHardProcess_, "genParticle_isHardProcess[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopy", &genParticle_isLastCopy_, "genParticle_isLastCopy[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopyBeforeFSR", &genParticle_isLastCopyBeforeFSR_, "genParticle_isLastCopyBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptDecayed", &genParticle_isPromptDecayed_, "genParticle_isPromptDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDecayedLeptonHadron", &genParticle_isDecayedLeptonHadron_, "genParticle_isDecayedLeptonHadron[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessBeforeFSR", &genParticle_fromHardProcessBeforeFSR_, "genParticle_fromHardProcessBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessDecayed", &genParticle_fromHardProcessDecayed_, "genParticle_fromHardProcessDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessFinalState", &genParticle_fromHardProcessFinalState_, "genParticle_fromHardProcessFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isMostlyLikePythia6Status3", &genParticle_isMostlyLikePythia6Status3_, "genParticle_isMostlyLikePythia6Status3[nGenParticle]/I");
  ntuple_->Branch("genParticle_l1pt", &genParticle_l1pt_, "genParticle_l1pt[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1eta", &genParticle_l1eta_, "genParticle_l1eta[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1phi", &genParticle_l1phi_, "genParticle_l1phi[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1charge", &genParticle_l1charge_, "genParticle_l1charge[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1q", &genParticle_l1q_, "genParticle_l1q[nGenParticle]/I");
  ntuple_->Branch("genParticle_l1dr", &genParticle_l1dr_, "genParticle_l1dr[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1ptByQ", &genParticle_l1ptByQ_, "genParticle_l1ptByQ[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1etaByQ", &genParticle_l1etaByQ_, "genParticle_l1etaByQ[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1phiByQ", &genParticle_l1phiByQ_, "genParticle_l1phiByQ[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1chargeByQ", &genParticle_l1chargeByQ_, "genParticle_l1chargeByQ[nGenParticle]/D");
  ntuple_->Branch("genParticle_l1qByQ", &genParticle_l1qByQ_, "genParticle_l1qByQ[nGenParticle]/I");
  ntuple_->Branch("genParticle_l1drByQ", &genParticle_l1drByQ_, "genParticle_l1drByQ[nGenParticle]/D");

  ntuple_->Branch("vec_firedTrigger", &vec_firedTrigger_);
  ntuple_->Branch("vec_filterName", &vec_filterName_);
  ntuple_->Branch("vec_HLTObj_pt", &vec_HLTObj_pt_);
  ntuple_->Branch("vec_HLTObj_eta", &vec_HLTObj_eta_);
  ntuple_->Branch("vec_HLTObj_phi", &vec_HLTObj_phi_);

  ntuple_->Branch("vec_myFiredTrigger", &vec_myFiredTrigger_);
  ntuple_->Branch("vec_myFilterName", &vec_myFilterName_);
  ntuple_->Branch("vec_myHLTObj_pt", &vec_myHLTObj_pt_);
  ntuple_->Branch("vec_myHLTObj_eta", &vec_myHLTObj_eta_);
  ntuple_->Branch("vec_myHLTObj_phi", &vec_myHLTObj_phi_);

  ntuple_->Branch("nMuon", &nMuon_, "nMuon/I");

  ntuple_->Branch("muon_pt", &muon_pt_, "muon_pt[nMuon]/D");
  ntuple_->Branch("muon_eta", &muon_eta_, "muon_eta[nMuon]/D");
  ntuple_->Branch("muon_phi", &muon_phi_, "muon_phi[nMuon]/D");
  ntuple_->Branch("muon_px", &muon_px_, "muon_px[nMuon]/D");
  ntuple_->Branch("muon_py", &muon_py_, "muon_py[nMuon]/D");
  ntuple_->Branch("muon_pz", &muon_pz_, "muon_pz[nMuon]/D");
  ntuple_->Branch("muon_dB", &muon_dB_, "muon_dB[nMuon]/D");
  ntuple_->Branch("muon_charge", &muon_charge_, "muon_charge[nMuon]/D");
  ntuple_->Branch("muon_isGLB", &muon_isGLB_, "muon_isGLB[nMuon]/I");
  ntuple_->Branch("muon_isSTA", &muon_isSTA_, "muon_isSTA[nMuon]/I");
  ntuple_->Branch("muon_isTRK", &muon_isTRK_, "muon_isTRK[nMuon]/I");
  ntuple_->Branch("muon_isPF", &muon_isPF_, "muon_isPF[nMuon]/I");
  ntuple_->Branch("muon_isTight", &muon_isTight_, "muon_isTight[nMuon]/I");
  ntuple_->Branch("muon_isMedium", &muon_isMedium_, "muon_isMedium[nMuon]/I");
  ntuple_->Branch("muon_isLoose", &muon_isLoose_, "muon_isLoose[nMuon]/I");
  ntuple_->Branch("muon_isHighPt", &muon_isHighPt_, "muon_isHighPt[nMuon]/I");
  ntuple_->Branch("muon_isHighPtNew", &muon_isHighPtNew_, "muon_isHighPtNew[nMuon]/I");
  ntuple_->Branch("muon_isSoft", &muon_isSoft_, "muon_isSoft[nMuon]/I");

  ntuple_->Branch("muon_iso03_sumPt", &muon_iso03_sumPt_, "muon_iso03_sumPt[nMuon]/D");
  ntuple_->Branch("muon_iso03_hadEt", &muon_iso03_hadEt_, "muon_iso03_hadEt[nMuon]/D");
  ntuple_->Branch("muon_iso03_emEt", &muon_iso03_emEt_, "muon_iso03_emEt[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_charged", &muon_PFIso03_charged_, "muon_PFIso03_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_neutral", &muon_PFIso03_neutral_, "muon_PFIso03_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_photon", &muon_PFIso03_photon_, "muon_PFIso03_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_sumPU", &muon_PFIso03_sumPU_, "muon_PFIso03_sumPU[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_charged", &muon_PFIso04_charged_, "muon_PFIso04_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_neutral", &muon_PFIso04_neutral_, "muon_PFIso04_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_photon", &muon_PFIso04_photon_, "muon_PFIso04_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_sumPU", &muon_PFIso04_sumPU_, "muon_PFIso04_sumPU[nMuon]/D");

  ntuple_->Branch("muon_PFCluster03_ECAL", &muon_PFCluster03_ECAL_, "muon_PFCluster03_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster03_HCAL", &muon_PFCluster03_HCAL_, "muon_PFCluster03_HCAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_ECAL", &muon_PFCluster04_ECAL_, "muon_PFCluster04_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_HCAL", &muon_PFCluster04_HCAL_, "muon_PFCluster04_HCAL[nMuon]/D");
  ntuple_->Branch("muon_normChi2_global", &muon_normChi2_global_, "muon_normChi2_global[nMuon]/D");
  ntuple_->Branch("muon_nTrackerHit_global", &muon_nTrackerHit_global_, "muon_nTrackerHit_global[nMuon]/I");
  ntuple_->Branch("muon_nTrackerLayer_global", &muon_nTrackerLayer_global_, "muon_nTrackerLayer_global[nMuon]/I");
  ntuple_->Branch("muon_nPixelHit_global", &muon_nPixelHit_global_, "muon_nPixelHit_global[nMuon]/I");
  ntuple_->Branch("muon_nMuonHit_global", &muon_nMuonHit_global_, "muon_nMuonHit_global[nMuon]/I");
  ntuple_->Branch("muon_normChi2_inner", &muon_normChi2_inner_, "muon_normChi2_inner[nMuon]/D");
  ntuple_->Branch("muon_nTrackerHit_inner", &muon_nTrackerHit_inner_, "muon_nTrackerHit_inner[nMuon]/I");
  ntuple_->Branch("muon_nTrackerLayer_inner", &muon_nTrackerLayer_inner_, "muon_nTrackerLayer_inner[nMuon]/I");
  ntuple_->Branch("muon_nPixelHit_inner", &muon_nPixelHit_inner_, "muon_nPixelHit_inner[nMuon]/I");
  ntuple_->Branch("muon_pt_tuneP", &muon_pt_tuneP_, "muon_pt_tuneP[nMuon]/D");
  ntuple_->Branch("muon_ptError_tuneP", &muon_ptError_tuneP_, "muon_ptError_tuneP[nMuon]/D");
  ntuple_->Branch("muon_dxyVTX_best", &muon_dxyVTX_best_, "muon_dxyVTX_best[nMuon]/D");
  ntuple_->Branch("muon_dzVTX_best", &muon_dzVTX_best_, "muon_dzVTX_best[nMuon]/D");
  ntuple_->Branch("muon_nMatchedStation", &muon_nMatchedStation_, "muon_nMatchedStation[nMuon]/I");
  ntuple_->Branch("muon_nMatchedRPCLayer", &muon_nMatchedRPCLayer_, "muon_nMatchedRPCLayer[nMuon]/I");
  ntuple_->Branch("muon_stationMask", &muon_stationMask_, "muon_stationMask[nMuon]/I");
  ntuple_->Branch("muon_dxy_bs", &muon_dxy_bs_, "muon_dxy_bs[nMuon]/D");
  ntuple_->Branch("muon_dxyError_bs", &muon_dxyError_bs_, "muon_dxyError_bs[nMuon]/D");
  ntuple_->Branch("muon_dz_bs", &muon_dz_bs_, "muon_dz_bs[nMuon]/D");
  ntuple_->Branch("muon_dzError", &muon_dzError_, "muon_dzError[nMuon]/D");
  ntuple_->Branch("muon_IPSig", &muon_IPSig_, "muon_IPSig[nMuon]/D");

  ntuple_->Branch("nL3Muon", &nL3Muon_, "nL3Muon/I");
  ntuple_->Branch("L3Muon_pt", &L3Muon_pt_, "L3Muon_pt[nL3Muon]/D");
  ntuple_->Branch("L3Muon_eta", &L3Muon_eta_, "L3Muon_eta[nL3Muon]/D");
  ntuple_->Branch("L3Muon_phi", &L3Muon_phi_, "L3Muon_phi[nL3Muon]/D");
  ntuple_->Branch("L3Muon_charge", &L3Muon_charge_, "L3Muon_charge[nL3Muon]/D");
  ntuple_->Branch("L3Muon_trkPt", &L3Muon_trkPt_, "L3Muon_trkPt[nL3Muon]/D");

  ntuple_->Branch("nL2Muon", &nL2Muon_, "nL2Muon/I");
  ntuple_->Branch("L2Muon_pt", &L2Muon_pt_, "L2Muon_pt[nL2Muon]/D");
  ntuple_->Branch("L2Muon_eta", &L2Muon_eta_, "L2Muon_eta[nL2Muon]/D");
  ntuple_->Branch("L2Muon_phi", &L2Muon_phi_, "L2Muon_phi[nL2Muon]/D");
  ntuple_->Branch("L2Muon_charge", &L2Muon_charge_, "L2Muon_charge[nL2Muon]/D");
  ntuple_->Branch("L2Muon_trkPt", &L2Muon_trkPt_, "L2Muon_trkPt[nL2Muon]/D");

  ntuple_->Branch("nTkMuon", &nTkMuon_, "nTkMuon/I");
  ntuple_->Branch("TkMuon_pt", &TkMuon_pt_, "TkMuon_pt[nTkMuon]/D");
  ntuple_->Branch("TkMuon_eta", &TkMuon_eta_, "TkMuon_eta[nTkMuon]/D");
  ntuple_->Branch("TkMuon_phi", &TkMuon_phi_, "TkMuon_phi[nTkMuon]/D");
  ntuple_->Branch("TkMuon_charge", &TkMuon_charge_, "TkMuon_charge[nTkMuon]/D");
  ntuple_->Branch("TkMuon_trkPt", &TkMuon_trkPt_, "TkMuon_trkPt[nTkMuon]/D");

  ntuple_->Branch("nL1Muon", &nL1Muon_, "nL1Muon/I");
  ntuple_->Branch("L1Muon_pt", &L1Muon_pt_, "L1Muon_pt[nL1Muon]/D");
  ntuple_->Branch("L1Muon_eta", &L1Muon_eta_, "L1Muon_eta[nL1Muon]/D");
  ntuple_->Branch("L1Muon_phi", &L1Muon_phi_, "L1Muon_phi[nL1Muon]/D");
  ntuple_->Branch("L1Muon_charge", &L1Muon_charge_, "L1Muon_charge[nL1Muon]/D");
  ntuple_->Branch("L1Muon_quality", &L1Muon_quality_, "L1Muon_quality[nL1Muon]/D");
  ntuple_->Branch("L1Muon_etaAtVtx", &L1Muon_etaAtVtx_, "L1Muon_etaAtVtx[nL1Muon]/D");
  ntuple_->Branch("L1Muon_phiAtVtx", &L1Muon_phiAtVtx_, "L1Muon_phiAtVtx[nL1Muon]/D");

  ntuple_->Branch("nIterL3OI", &nIterL3OI_, "nIterL3OI/I");
  ntuple_->Branch("iterL3OI_inner_pt", &iterL3OI_inner_pt_, "iterL3OI_inner_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_eta", &iterL3OI_inner_eta_, "iterL3OI_inner_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_phi", &iterL3OI_inner_phi_, "iterL3OI_inner_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_charge", &iterL3OI_inner_charge_, "iterL3OI_inner_charge[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_pt", &iterL3OI_outer_pt_, "iterL3OI_outer_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_eta", &iterL3OI_outer_eta_, "iterL3OI_outer_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_phi", &iterL3OI_outer_phi_, "iterL3OI_outer_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_charge", &iterL3OI_outer_charge_, "iterL3OI_outer_charge[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_pt", &iterL3OI_global_pt_, "iterL3OI_global_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_eta", &iterL3OI_global_eta_, "iterL3OI_global_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_phi", &iterL3OI_global_phi_, "iterL3OI_global_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_charge", &iterL3OI_global_charge_, "iterL3OI_global_charge[nIterL3OI]/D");

  ntuple_->Branch("nIterL3IOFromL2", &nIterL3IOFromL2_, "nIterL3IOFromL2/I");
  ntuple_->Branch("iterL3IOFromL2_inner_pt", &iterL3IOFromL2_inner_pt_, "iterL3IOFromL2_inner_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_eta", &iterL3IOFromL2_inner_eta_, "iterL3IOFromL2_inner_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_phi", &iterL3IOFromL2_inner_phi_, "iterL3IOFromL2_inner_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_charge", &iterL3IOFromL2_inner_charge_, "iterL3IOFromL2_inner_charge[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_pt", &iterL3IOFromL2_outer_pt_, "iterL3IOFromL2_outer_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_eta", &iterL3IOFromL2_outer_eta_, "iterL3IOFromL2_outer_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_phi", &iterL3IOFromL2_outer_phi_, "iterL3IOFromL2_outer_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_charge", &iterL3IOFromL2_outer_charge_, "iterL3IOFromL2_outer_charge[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_pt", &iterL3IOFromL2_global_pt_, "iterL3IOFromL2_global_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_eta", &iterL3IOFromL2_global_eta_, "iterL3IOFromL2_global_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_phi", &iterL3IOFromL2_global_phi_, "iterL3IOFromL2_global_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_charge", &iterL3IOFromL2_global_charge_, "iterL3IOFromL2_global_charge[nIterL3IOFromL2]/D");

  ntuple_->Branch("nIterL3IOFromL1", &nIterL3IOFromL1_, "nIterL3IOFromL1/I");
  ntuple_->Branch("iterL3IOFromL1_pt", &iterL3IOFromL1_pt_, "iterL3IOFromL1_pt[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_eta", &iterL3IOFromL1_eta_, "iterL3IOFromL1_eta[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_phi", &iterL3IOFromL1_phi_, "iterL3IOFromL1_phi[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_charge", &iterL3IOFromL1_charge_, "iterL3IOFromL1_charge[nIterL3IOFromL1]/D");

  ntuple_->Branch("nIterL3FromL2", &nIterL3FromL2_, "nIterL3FromL2/I");
  ntuple_->Branch("iterL3FromL2_inner_pt", &iterL3FromL2_inner_pt_, "iterL3FromL2_inner_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_eta", &iterL3FromL2_inner_eta_, "iterL3FromL2_inner_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_phi", &iterL3FromL2_inner_phi_, "iterL3FromL2_inner_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_charge", &iterL3FromL2_inner_charge_, "iterL3FromL2_inner_charge[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_pt", &iterL3FromL2_outer_pt_, "iterL3FromL2_outer_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_eta", &iterL3FromL2_outer_eta_, "iterL3FromL2_outer_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_phi", &iterL3FromL2_outer_phi_, "iterL3FromL2_outer_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_charge", &iterL3FromL2_outer_charge_, "iterL3FromL2_outer_charge[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_pt", &iterL3FromL2_global_pt_, "iterL3FromL2_global_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_eta", &iterL3FromL2_global_eta_, "iterL3FromL2_global_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_phi", &iterL3FromL2_global_phi_, "iterL3FromL2_global_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_charge", &iterL3FromL2_global_charge_, "iterL3FromL2_global_charge[nIterL3FromL2]/D");

  ntuple_->Branch("nIterL3MuonNoID",       &nIterL3MuonNoID_,       "nIterL3MuonNoID/I");
  ntuple_->Branch("iterL3MuonNoID_pt",     &iterL3MuonNoID_pt_,     "iterL3MuonNoID_pt[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_innerPt",     &iterL3MuonNoID_innerPt_,     "iterL3MuonNoID_innerPt[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_eta",    &iterL3MuonNoID_eta_,    "iterL3MuonNoID_eta[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_phi",    &iterL3MuonNoID_phi_,    "iterL3MuonNoID_phi[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_charge", &iterL3MuonNoID_charge_, "iterL3MuonNoID_charge[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_isGLB",  &iterL3MuonNoID_isGLB_,  "iterL3MuonNoID_isGLB[nIterL3MuonNoID]/I");
  ntuple_->Branch("iterL3MuonNoID_isSTA",  &iterL3MuonNoID_isSTA_,  "iterL3MuonNoID_isSTA[nIterL3MuonNoID]/I");
  ntuple_->Branch("iterL3MuonNoID_isTRK",  &iterL3MuonNoID_isTRK_,  "iterL3MuonNoID_isTRK[nIterL3MuonNoID]/I");

  ntuple_->Branch("nIterL3Muon",       &nIterL3Muon_,       "nIterL3Muon/I");
  ntuple_->Branch("iterL3Muon_pt",     &iterL3Muon_pt_,     "iterL3Muon_pt[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_innerPt", &iterL3Muon_innerPt_, "iterL3Muon_innerPt[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_eta",    &iterL3Muon_eta_,    "iterL3Muon_eta[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_phi",    &iterL3Muon_phi_,    "iterL3Muon_phi[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_charge", &iterL3Muon_charge_, "iterL3Muon_charge[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_isGLB",  &iterL3Muon_isGLB_,  "iterL3Muon_isGLB[nIterL3Muon]/I");
  ntuple_->Branch("iterL3Muon_isSTA",  &iterL3Muon_isSTA_,  "iterL3Muon_isSTA[nIterL3Muon]/I");
  ntuple_->Branch("iterL3Muon_isTRK",  &iterL3Muon_isTRK_,  "iterL3Muon_isTRK[nIterL3Muon]/I");

  //SThltIterL3OISeedsFromL2Muons->setBranch(ntuple_,"hltIterL3OISeedsFromL2Muons");
  //SThltIter0IterL3MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3MuonPixelSeedsFromPixelTracks");
  //SThltIter2IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3MuonPixelSeeds");
  //SThltIter3IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3MuonPixelSeeds");
  //SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks");
  //SThltIter2IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3FromL1MuonPixelSeeds");
  //SThltIter3IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3FromL1MuonPixelSeeds");

  TThltIterL3OIMuonTrack->setBranch(ntuple_,"hltIterL3OIMuonTrack");
  TThltIter0IterL3MuonTrack->setBranch(ntuple_,"hltIter0IterL3MuonTrack");
  TThltIter2IterL3MuonTrack->setBranch(ntuple_,"hltIter2IterL3MuonTrack");
  TThltIter3IterL3MuonTrack->setBranch(ntuple_,"hltIter3IterL3MuonTrack");
  TThltIter0IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter0IterL3FromL1MuonTrack");
  TThltIter2IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter2IterL3FromL1MuonTrack");
  TThltIter3IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter3IterL3FromL1MuonTrack");

  TrkParticle->setBranch(ntuple_,"TP");

  VThltIterL3MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3MuonTrimmedPixelVertices");
  VThltIterL3FromL1MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3FromL1MuonTrimmedPixelVertices");

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    TString trkName = TString(trackCollectionNames_.at(i));
    TString tpName  = "tpTo_" + TString(trackCollectionNames_.at(i));

    trkTemplates_.at(i)->setBranch(ntuple_, trkName );
    tpTemplates_.at(i)->setBranch(ntuple_,  tpName );
  }

}

void MuonHLTNtupler::Fill_Muon(const edm::Event &iEvent)
{
  edm::Handle<std::vector<reco::Muon> > h_offlineMuon;
  if( iEvent.getByToken(t_offlineMuon_, h_offlineMuon) ) // -- only when the dataset has offline muon collection (e.g. AOD) -- //
  {
    edm::Handle<reco::VertexCollection> h_offlineVertex;
    iEvent.getByToken(t_offlineVertex_, h_offlineVertex);
    const reco::Vertex & pv = h_offlineVertex->at(0);

    int _nMuon = 0;
    for(std::vector<reco::Muon>::const_iterator mu=h_offlineMuon->begin(); mu!=h_offlineMuon->end(); ++mu)
    {
      muon_pt_[_nMuon]  = mu->pt();
      muon_eta_[_nMuon] = mu->eta();
      muon_phi_[_nMuon] = mu->phi();
      muon_px_[_nMuon]  = mu->px();
      muon_py_[_nMuon]  = mu->py();
      muon_pz_[_nMuon]  = mu->pz();
      muon_charge_[_nMuon] = mu->charge();

      if( mu->isGlobalMuon() ) muon_isGLB_[_nMuon] = 1;
      if( mu->isStandAloneMuon() ) muon_isSTA_[_nMuon] = 1;
      if( mu->isTrackerMuon() ) muon_isTRK_[_nMuon] = 1;
      if( mu->isPFMuon() ) muon_isPF_[_nMuon] = 1;

      // -- defintion of ID functions: http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_0/doc/html/da/d18/namespacemuon.html#ac122b2516e5711ce206256d7945473d2 -- //
      if( muon::isTightMuon( (*mu), pv ) )  muon_isTight_[_nMuon] = 1;
      if( muon::isMediumMuon( (*mu) ) )     muon_isMedium_[_nMuon] = 1;
      if( muon::isLooseMuon( (*mu) ) )      muon_isLoose_[_nMuon] = 1;
      if( muon::isHighPtMuon( (*mu), pv ) ) muon_isHighPt_[_nMuon] = 1;
      if( isNewHighPtMuon( (*mu), pv ) )    muon_isHighPtNew_[_nMuon] = 1;

      // -- bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx, bool run2016_hip_mitigation)
      // -- it is different under CMSSW_8_0_29: bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx)
      // -- Remove this part to avoid compile error (and soft muon would not be used for now) - need to be fixed at some point
      // if( muon::isSoftMuon( (*mu), pv, 0) ) muon_isSoft_[_nMuon] = 1;

      muon_iso03_sumPt_[_nMuon] = mu->isolationR03().sumPt;
      muon_iso03_hadEt_[_nMuon] = mu->isolationR03().hadEt;
      muon_iso03_emEt_[_nMuon]  = mu->isolationR03().emEt;

      muon_PFIso03_charged_[_nMuon] = mu->pfIsolationR03().sumChargedHadronPt;
      muon_PFIso03_neutral_[_nMuon] = mu->pfIsolationR03().sumNeutralHadronEt;
      muon_PFIso03_photon_[_nMuon]  = mu->pfIsolationR03().sumPhotonEt;
      muon_PFIso03_sumPU_[_nMuon]   = mu->pfIsolationR03().sumPUPt;

      muon_PFIso04_charged_[_nMuon] = mu->pfIsolationR04().sumChargedHadronPt;
      muon_PFIso04_neutral_[_nMuon] = mu->pfIsolationR04().sumNeutralHadronEt;
      muon_PFIso04_photon_[_nMuon]  = mu->pfIsolationR04().sumPhotonEt;
      muon_PFIso04_sumPU_[_nMuon]   = mu->pfIsolationR04().sumPUPt;

      // reco::MuonRef muRef = reco::MuonRef(h_offlineMuon, _nMuon);

      reco::TrackRef globalTrk = mu->globalTrack();
      if( globalTrk.isNonnull() )
      {
        muon_normChi2_global_[_nMuon] = globalTrk->normalizedChi2();

        const reco::HitPattern & globalTrkHit = globalTrk->hitPattern();
        muon_nTrackerHit_global_[_nMuon]   = globalTrkHit.numberOfValidTrackerHits();
        muon_nTrackerLayer_global_[_nMuon] = globalTrkHit.trackerLayersWithMeasurement();
        muon_nPixelHit_global_[_nMuon]     = globalTrkHit.numberOfValidPixelHits();
        muon_nMuonHit_global_[_nMuon]      = globalTrkHit.numberOfValidMuonHits();
      }

      reco::TrackRef innerTrk = mu->innerTrack();
      if( innerTrk.isNonnull() )
      {
        muon_normChi2_inner_[_nMuon] = innerTrk->normalizedChi2();

        muon_dxy_bs_[_nMuon] = innerTrk->dxy(bs->position());
        muon_dxyError_bs_[_nMuon] = innerTrk->dxyError(*bs);
        muon_dz_bs_[_nMuon] = innerTrk->dz(bs->position());
        muon_dzError_[_nMuon] = innerTrk->dzError();
        if (innerTrk->dxyError(*bs) > 0.) {
          muon_IPSig_[_nMuon] = abs(innerTrk->dxy(bs->position()) / innerTrk->dxyError(*bs));
        }

        const reco::HitPattern & innerTrkHit = innerTrk->hitPattern();
        muon_nTrackerHit_inner_[_nMuon]   = innerTrkHit.numberOfValidTrackerHits();
        muon_nTrackerLayer_inner_[_nMuon] = innerTrkHit.trackerLayersWithMeasurement();
        muon_nPixelHit_inner_[_nMuon]     = innerTrkHit.numberOfValidPixelHits();
      }

      reco::TrackRef tunePTrk = mu->tunePMuonBestTrack();
      if( tunePTrk.isNonnull() )
      {
        muon_pt_tuneP_[_nMuon]      = tunePTrk->pt();
        muon_ptError_tuneP_[_nMuon] = tunePTrk->ptError();
      }

      muon_dxyVTX_best_[_nMuon] = mu->muonBestTrack()->dxy( pv.position() );
      muon_dzVTX_best_[_nMuon]  = mu->muonBestTrack()->dz( pv.position() );

      muon_nMatchedStation_[_nMuon] = mu->numberOfMatchedStations();
      muon_nMatchedRPCLayer_[_nMuon] = mu->numberOfMatchedRPCLayers();
      muon_stationMask_[_nMuon] = mu->stationMask();

      _nMuon++;
    }

    nMuon_ = _nMuon;
  }
}

void MuonHLTNtupler::Fill_HLT(const edm::Event &iEvent, bool isMYHLT)
{
  edm::Handle<edm::TriggerResults>  h_triggerResults;
  edm::Handle<trigger::TriggerEvent> h_triggerEvent;

  if( isMYHLT )
  {
    iEvent.getByToken(t_myTriggerResults_, h_triggerResults);
    iEvent.getByToken(t_myTriggerEvent_,   h_triggerEvent);
  }
  else
  {
    iEvent.getByToken(t_triggerResults_, h_triggerResults);
    iEvent.getByToken(t_triggerEvent_,   h_triggerEvent);
  }

  edm::TriggerNames triggerNames = iEvent.triggerNames(*h_triggerResults);

  for(unsigned int itrig=0; itrig<triggerNames.size(); ++itrig)
  {
    LogDebug("triggers") << triggerNames.triggerName(itrig);

    if( h_triggerResults->accept(itrig) )
    {
      std::string pathName = triggerNames.triggerName(itrig);
      if( SavedTriggerCondition(pathName) )
      {
        if( isMYHLT ) vec_myFiredTrigger_.push_back( pathName );
        else          vec_firedTrigger_.push_back( pathName );
      }
    } // -- end of if fired -- //

  } // -- end of iteration over all trigger names -- //

  const trigger::size_type nFilter(h_triggerEvent->sizeFilters());
  for( trigger::size_type i_filter=0; i_filter<nFilter; i_filter++)
  {
    std::string filterName = h_triggerEvent->filterTag(i_filter).encode();

    if( SavedFilterCondition(filterName) )
    {
      trigger::Keys objectKeys = h_triggerEvent->filterKeys(i_filter);
      const trigger::TriggerObjectCollection& triggerObjects(h_triggerEvent->getObjects());

      for( trigger::size_type i_key=0; i_key<objectKeys.size(); i_key++)
      {
        trigger::size_type objKey = objectKeys.at(i_key);
        const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

        if( isMYHLT )
        {
          vec_myFilterName_.push_back( filterName );
          vec_myHLTObj_pt_.push_back( triggerObj.pt() );
          vec_myHLTObj_eta_.push_back( triggerObj.eta() );
          vec_myHLTObj_phi_.push_back( triggerObj.phi() );
        }
        else
        {
          vec_filterName_.push_back( filterName );
          vec_HLTObj_pt_.push_back( triggerObj.pt() );
          vec_HLTObj_eta_.push_back( triggerObj.eta() );
          vec_HLTObj_phi_.push_back( triggerObj.phi() );
        }
      }
    } // -- end of if( muon filters )-- //
  } // -- end of filter iteration -- //
}

bool MuonHLTNtupler::SavedTriggerCondition( std::string& pathName )
{
  bool flag = false;

  // -- muon triggers
  if( pathName.find("Mu")           != std::string::npos ||
      pathName.find("HLT_IsoMu")    != std::string::npos ||
      pathName.find("HLT_Mu")       != std::string::npos ||
      pathName.find("HLT_OldMu")    != std::string::npos ||
      pathName.find("HLT_TkMu")     != std::string::npos ||
      pathName.find("HLT_IsoTkMu")  != std::string::npos ||
      pathName.find("HLT_DoubleMu") != std::string::npos ||
      pathName.find("HLT_Mu8_T")    != std::string::npos ) flag = true;

  return flag;
}

bool MuonHLTNtupler::SavedFilterCondition( std::string& filterName )
{
  bool flag = false;

  // -- muon filters
  if( (
        filterName.find("sMu") != std::string::npos ||
        filterName.find("SingleMu") != std::string::npos ||
        filterName.find("TkMu") != std::string::npos
      ) &&
       filterName.find("Tau")      == std::string::npos &&
       filterName.find("EG")       == std::string::npos &&
       filterName.find("MultiFit") == std::string::npos ) flag = true;

  return flag;
}

void MuonHLTNtupler::Fill_HLTMuon(const edm::Event &iEvent)
{
  ///////////////////
  // -- L3 Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_L3Muon;
  if( iEvent.getByToken( t_L3Muon_, h_L3Muon ) )
  {
    int _nL3Muon = 0;
    for(unsigned int i_L3=0; i_L3<h_L3Muon->size(); i_L3++)
    {
      reco::RecoChargedCandidateRef ref_L3Mu(h_L3Muon, _nL3Muon);

      L3Muon_pt_[_nL3Muon]     = ref_L3Mu->pt();
      L3Muon_eta_[_nL3Muon]    = ref_L3Mu->eta();
      L3Muon_phi_[_nL3Muon]    = ref_L3Mu->phi();
      L3Muon_charge_[_nL3Muon] = ref_L3Mu->charge();

      reco::TrackRef trackRef = ref_L3Mu->track();
      L3Muon_trkPt_[_nL3Muon] = trackRef->pt();

      _nL3Muon++;
    }
    nL3Muon_ = _nL3Muon;
  } // -- if( L3 handle is valid ) -- //


  ///////////////////
  // -- L2 Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  if( iEvent.getByToken( t_L2Muon_, h_L2Muon ) )
  {
    int _nL2Muon = 0;
    for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
    {
      reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, _nL2Muon);

      L2Muon_pt_[_nL2Muon]     = ref_L2Mu->pt();
      L2Muon_eta_[_nL2Muon]    = ref_L2Mu->eta();
      L2Muon_phi_[_nL2Muon]    = ref_L2Mu->phi();
      L2Muon_charge_[_nL2Muon] = ref_L2Mu->charge();

      reco::TrackRef trackRef = ref_L2Mu->track();
      L2Muon_trkPt_[_nL2Muon] = trackRef->pt();

      _nL2Muon++;
    }
    nL2Muon_ = _nL2Muon;
  }

  ///////////////////
  // -- Tk Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_TkMuon;
  if( iEvent.getByToken( t_TkMuon_, h_TkMuon ) )
  {
    int _nTkMuon = 0;
    for( unsigned int i_Tk=0; i_Tk<h_TkMuon->size(); i_Tk++)
    {
      reco::RecoChargedCandidateRef ref_TkMu(h_TkMuon, _nTkMuon);

      TkMuon_pt_[_nTkMuon]     = ref_TkMu->pt();
      TkMuon_eta_[_nTkMuon]    = ref_TkMu->eta();
      TkMuon_phi_[_nTkMuon]    = ref_TkMu->phi();
      TkMuon_charge_[_nTkMuon] = ref_TkMu->charge();

      reco::TrackRef trackRef = ref_TkMu->track();
      TkMuon_trkPt_[_nTkMuon] = trackRef->pt();

      _nTkMuon++;
    }
    nTkMuon_ = _nTkMuon;
  }
}

void MuonHLTNtupler::Fill_L1Muon(const edm::Event &iEvent)
{
  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  if( iEvent.getByToken(t_L1Muon_, h_L1Muon) )
  {
    int _nL1Muon = 0;
    for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
    {
      if(ibx != 0) continue; // -- only take when ibx == 0 -- //
      for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
      {
        l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

        L1Muon_pt_[_nL1Muon]      = ref_L1Mu->pt();
        L1Muon_eta_[_nL1Muon]     = ref_L1Mu->eta();
        L1Muon_phi_[_nL1Muon]     = ref_L1Mu->phi();
        L1Muon_charge_[_nL1Muon]  = ref_L1Mu->charge();
        L1Muon_quality_[_nL1Muon] = ref_L1Mu->hwQual();

        L1Muon_etaAtVtx_[_nL1Muon] = ref_L1Mu->etaAtVtx();
        L1Muon_phiAtVtx_[_nL1Muon] = ref_L1Mu->phiAtVtx();

        _nL1Muon++;
      }
    }
    nL1Muon_ = _nL1Muon;
  }
}

void MuonHLTNtupler::Fill_GenParticle(const edm::Event &iEvent)
{
  // -- Gen-weight info -- //
  edm::Handle<GenEventInfoProduct> h_genEventInfo;
  iEvent.getByToken(t_genEventInfo_, h_genEventInfo);
  genEventWeight_ = h_genEventInfo->weight();

  // -- Gen-particle info -- //
  edm::Handle<edm::View<reco::GenParticle>> h_genParticle;
  iEvent.getByToken(t_genParticle_, h_genParticle);
  edm::Handle<pat::TriggerObjectStandAloneMatch> h_l1Matches;
  iEvent.getByToken(t_l1Matches_, h_l1Matches);
  edm::Handle<edm::ValueMap<int>> h_l1Qualities;
  iEvent.getByToken(t_l1MatchesQuality_, h_l1Qualities);
  edm::Handle<edm::ValueMap<float>> h_l1Drs;
  iEvent.getByToken(t_l1MatchesDeltaR_, h_l1Drs);
  edm::Handle<pat::TriggerObjectStandAloneMatch> h_l1MatchesByQ;
  iEvent.getByToken(t_l1MatchesByQ_, h_l1MatchesByQ);
  edm::Handle<edm::ValueMap<int>> h_l1QualitiesByQ;
  iEvent.getByToken(t_l1MatchesByQQuality_, h_l1QualitiesByQ);
  edm::Handle<edm::ValueMap<float>> h_l1DrsByQ;
  iEvent.getByToken(t_l1MatchesByQDeltaR_, h_l1DrsByQ);

  int _nGenParticle = 0;
  for( size_t i=0; i< h_genParticle->size(); ++i)
  {
    const auto &parCand = (*h_genParticle)[i];
    auto genRef = h_genParticle->refAt(i);

    if( abs(parCand.pdgId()) == 13 ) // -- only muons -- //
    {
      genParticle_ID_[_nGenParticle]     = parCand.pdgId();
      genParticle_status_[_nGenParticle] = parCand.status();
      genParticle_mother_[_nGenParticle] = parCand.mother(0)? parCand.mother(0)->pdgId(): -999;

      genParticle_pt_[_nGenParticle]  = parCand.pt();
      genParticle_eta_[_nGenParticle] = parCand.eta();
      genParticle_phi_[_nGenParticle] = parCand.phi();
      genParticle_px_[_nGenParticle]  = parCand.px();
      genParticle_py_[_nGenParticle]  = parCand.py();
      genParticle_pz_[_nGenParticle]  = parCand.pz();
      genParticle_energy_[_nGenParticle] = parCand.energy();
      genParticle_charge_[_nGenParticle] = parCand.charge();

      if( parCand.statusFlags().isPrompt() )                genParticle_isPrompt_[_nGenParticle] = 1;
      if( parCand.statusFlags().isTauDecayProduct() )       genParticle_isTauDecayProduct_[_nGenParticle] = 1;
      if( parCand.statusFlags().isPromptTauDecayProduct() ) genParticle_isPromptTauDecayProduct_[_nGenParticle] = 1;
      if( parCand.statusFlags().isDecayedLeptonHadron() )   genParticle_isDecayedLeptonHadron_[_nGenParticle] = 1;

      if( parCand.isPromptFinalState() ) genParticle_isPromptFinalState_[_nGenParticle] = 1;
      if( parCand.isDirectPromptTauDecayProductFinalState() ) genParticle_isDirectPromptTauDecayProductFinalState_[_nGenParticle] = 1;
      if( parCand.isHardProcess() ) genParticle_isHardProcess_[_nGenParticle] = 1;
      if( parCand.isLastCopy() ) genParticle_isLastCopy_[_nGenParticle] = 1;
      if( parCand.isLastCopyBeforeFSR() ) genParticle_isLastCopyBeforeFSR_[_nGenParticle] = 1;

      if( parCand.isPromptDecayed() )           genParticle_isPromptDecayed_[_nGenParticle] = 1;
      if( parCand.fromHardProcessBeforeFSR() )  genParticle_fromHardProcessBeforeFSR_[_nGenParticle] = 1;
      if( parCand.fromHardProcessDecayed() )    genParticle_fromHardProcessDecayed_[_nGenParticle] = 1;
      if( parCand.fromHardProcessFinalState() ) genParticle_fromHardProcessFinalState_[_nGenParticle] = 1;
      // if( parCand.isMostlyLikePythia6Status3() ) this->genParticle_isMostlyLikePythia6Status3[_nGenParticle] = 1;

      pat::TriggerObjectStandAloneRef genl1Match = (*h_l1Matches)[genRef];
      if (genl1Match.isNonnull()) {
        genParticle_l1pt_[_nGenParticle]      = genl1Match->pt();
        genParticle_l1eta_[_nGenParticle]     = genl1Match->eta();
        genParticle_l1phi_[_nGenParticle]     = genl1Match->phi();
        genParticle_l1charge_[_nGenParticle]  = genl1Match->charge();
        genParticle_l1q_[_nGenParticle]       = (*h_l1Qualities)[genRef];
        genParticle_l1dr_[_nGenParticle]      = (*h_l1Drs)[genRef];
      }

      pat::TriggerObjectStandAloneRef genl1MatchByQ = (*h_l1MatchesByQ)[genRef];
      if (genl1MatchByQ.isNonnull()) {
        genParticle_l1ptByQ_[_nGenParticle]      = genl1MatchByQ->pt();
        genParticle_l1etaByQ_[_nGenParticle]     = genl1MatchByQ->eta();
        genParticle_l1phiByQ_[_nGenParticle]     = genl1MatchByQ->phi();
        genParticle_l1chargeByQ_[_nGenParticle]  = genl1MatchByQ->charge();
        genParticle_l1qByQ_[_nGenParticle]       = (*h_l1QualitiesByQ)[genRef];
        genParticle_l1drByQ_[_nGenParticle]      = (*h_l1DrsByQ)[genRef];
      }

      _nGenParticle++;
    }
  }
  nGenParticle_ = _nGenParticle;
}

void MuonHLTNtupler::Fill_IterL3(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //////////////////////////
  // -- IterL3Muon -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Muon> > h_iterL3Muon;
  if( iEvent.getByToken( t_iterL3Muon_, h_iterL3Muon) )
  {
    int _nIterL3Muon = 0;
    for( auto i=0U; i<h_iterL3Muon->size(); ++i )
    {
      const auto& muon(h_iterL3Muon->at(i));

      iterL3Muon_pt_[_nIterL3Muon]     = muon.pt();
      iterL3Muon_eta_[_nIterL3Muon]    = muon.eta();
      iterL3Muon_phi_[_nIterL3Muon]    = muon.phi();
      iterL3Muon_charge_[_nIterL3Muon] = muon.charge();

      if( muon.isGlobalMuon() )     iterL3Muon_isGLB_[_nIterL3Muon] = 1;
      if( muon.isStandAloneMuon() ) iterL3Muon_isSTA_[_nIterL3Muon] = 1;
      if( muon.isTrackerMuon() )    iterL3Muon_isTRK_[_nIterL3Muon] = 1;

      reco::TrackRef innerTrk = muon.innerTrack();
      if( innerTrk.isNonnull() ) {
        iterL3Muon_innerPt_[_nIterL3Muon] = innerTrk->pt();

        tmpTrk trkTmp(innerTrk);
        iterL3IDpassed.push_back(trkTmp);

        const PTrajectoryStateOnDet tmpseed = innerTrk->seedRef()->startingState();
        tmpTSOD tsod(tmpseed);
        MuonIterSeedMap.insert(make_pair(tsod,i));
      }
      else {
        cout << "IterL3Muon: innerTrk.isNonnull(): this should never happen" << endl;
        tmpTrk trkTmp( -99999. );  // dummy tmpTrk
        iterL3IDpassed.push_back(trkTmp);
      }

      _nIterL3Muon++;
    } // -- end of muon iteration

    nIterL3Muon_ = _nIterL3Muon;
  } // -- if getByToken is valid

  ////////////////////
  // -- IterL3OI -- //
  ////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3OI;
  if( iEvent.getByToken( t_iterL3OI_, h_iterL3OI ) )
  {
    int _nIterL3OI = 0;
    for( unsigned int i=0; i<h_iterL3OI->size(); i++)
    {
      if( h_iterL3OI->at(i).trackerTrack().isNonnull() )
      {
        iterL3OI_inner_pt_[_nIterL3OI]     = h_iterL3OI->at(i).trackerTrack()->pt();
        iterL3OI_inner_eta_[_nIterL3OI]    = h_iterL3OI->at(i).trackerTrack()->eta();
        iterL3OI_inner_phi_[_nIterL3OI]    = h_iterL3OI->at(i).trackerTrack()->phi();
        iterL3OI_inner_charge_[_nIterL3OI] = h_iterL3OI->at(i).trackerTrack()->charge();
      }
      if( h_iterL3OI->at(i).standAloneTrack().isNonnull() )
      {
        iterL3OI_outer_pt_[_nIterL3OI]     = h_iterL3OI->at(i).standAloneTrack()->pt();
        iterL3OI_outer_eta_[_nIterL3OI]    = h_iterL3OI->at(i).standAloneTrack()->eta();
        iterL3OI_outer_phi_[_nIterL3OI]    = h_iterL3OI->at(i).standAloneTrack()->phi();
        iterL3OI_outer_charge_[_nIterL3OI] = h_iterL3OI->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3OI->at(i).globalTrack().isNonnull() )
      {
        iterL3OI_global_pt_[_nIterL3OI]     = h_iterL3OI->at(i).globalTrack()->pt();
        iterL3OI_global_eta_[_nIterL3OI]    = h_iterL3OI->at(i).globalTrack()->eta();
        iterL3OI_global_phi_[_nIterL3OI]    = h_iterL3OI->at(i).globalTrack()->phi();
        iterL3OI_global_charge_[_nIterL3OI] = h_iterL3OI->at(i).globalTrack()->charge();
      }
      _nIterL3OI++;
    }
    nIterL3OI_ = _nIterL3OI;
  }

  //////////////////////////
  // -- IterL3IOFromL2 -- //
  //////////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3IOFromL2;
  if( iEvent.getByToken( t_iterL3IOFromL2_, h_iterL3IOFromL2 ) )
  {
    int _nIterL3IOFromL2 = 0;
    for( unsigned int i=0; i<h_iterL3IOFromL2->size(); i++)
    {
      if( h_iterL3IOFromL2->at(i).trackerTrack().isNonnull() )
      {
        iterL3IOFromL2_inner_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).trackerTrack()->pt();
        iterL3IOFromL2_inner_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).trackerTrack()->eta();
        iterL3IOFromL2_inner_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).trackerTrack()->phi();
        iterL3IOFromL2_inner_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).trackerTrack()->charge();
      }
      if( h_iterL3IOFromL2->at(i).standAloneTrack().isNonnull() )
      {
        iterL3IOFromL2_outer_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).standAloneTrack()->pt();
        iterL3IOFromL2_outer_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).standAloneTrack()->eta();
        iterL3IOFromL2_outer_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).standAloneTrack()->phi();
        iterL3IOFromL2_outer_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3IOFromL2->at(i).globalTrack().isNonnull() )
      {
        iterL3IOFromL2_global_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).globalTrack()->pt();
        iterL3IOFromL2_global_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).globalTrack()->eta();
        iterL3IOFromL2_global_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).globalTrack()->phi();
        iterL3IOFromL2_global_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).globalTrack()->charge();
      }
      _nIterL3IOFromL2++;
    }
    nIterL3IOFromL2_ = _nIterL3IOFromL2;
  }

  ////////////////////////////////
  // -- IterL3FromL2 (OI+IO) -- //
  ////////////////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3FromL2;
  if( iEvent.getByToken( t_iterL3FromL2_, h_iterL3FromL2 ) )
  {
    int _nIterL3FromL2 = 0;
    for( unsigned int i=0; i<h_iterL3FromL2->size(); i++)
    {
      if( h_iterL3FromL2->at(i).trackerTrack().isNonnull() )
      {
        iterL3FromL2_inner_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).trackerTrack()->pt();
        iterL3FromL2_inner_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).trackerTrack()->eta();
        iterL3FromL2_inner_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).trackerTrack()->phi();
        iterL3FromL2_inner_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).trackerTrack()->charge();
      }
      if( h_iterL3FromL2->at(i).standAloneTrack().isNonnull() )
      {
        iterL3FromL2_outer_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).standAloneTrack()->pt();
        iterL3FromL2_outer_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).standAloneTrack()->eta();
        iterL3FromL2_outer_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).standAloneTrack()->phi();
        iterL3FromL2_outer_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3FromL2->at(i).globalTrack().isNonnull() )
      {
        iterL3FromL2_global_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).globalTrack()->pt();
        iterL3FromL2_global_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).globalTrack()->eta();
        iterL3FromL2_global_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).globalTrack()->phi();
        iterL3FromL2_global_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).globalTrack()->charge();
      }
      _nIterL3FromL2++;
    }
    nIterL3FromL2_ = _nIterL3FromL2;
  }

  //////////////////////////
  // -- IterL3IOFromL1 -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Track> > h_iterL3IOFromL1;
  if( iEvent.getByToken( t_iterL3IOFromL1_, h_iterL3IOFromL1 ) )
  {
    int _nIterL3IOFromL1 = 0;
    for( unsigned int i=0; i<h_iterL3IOFromL1->size(); i++)
    {
      iterL3IOFromL1_pt_[_nIterL3IOFromL1]     = h_iterL3IOFromL1->at(i).pt();
      iterL3IOFromL1_eta_[_nIterL3IOFromL1]    = h_iterL3IOFromL1->at(i).eta();
      iterL3IOFromL1_phi_[_nIterL3IOFromL1]    = h_iterL3IOFromL1->at(i).phi();
      iterL3IOFromL1_charge_[_nIterL3IOFromL1] = h_iterL3IOFromL1->at(i).charge();
      _nIterL3IOFromL1++;
    }
    nIterL3IOFromL1_ = _nIterL3IOFromL1;
  }

  //////////////////////////
  // -- IterL3MuonNoID -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Muon> > h_iterL3MuonNoID;
  if( iEvent.getByToken( t_iterL3MuonNoID_, h_iterL3MuonNoID) )
  {
    int _nIterL3MuonNoID = 0;
    for( auto i=0U; i<h_iterL3MuonNoID->size(); ++i )
    {
      const auto& muon(h_iterL3MuonNoID->at(i));

      iterL3MuonNoID_pt_[_nIterL3MuonNoID]     = muon.pt();
      iterL3MuonNoID_eta_[_nIterL3MuonNoID]    = muon.eta();
      iterL3MuonNoID_phi_[_nIterL3MuonNoID]    = muon.phi();
      iterL3MuonNoID_charge_[_nIterL3MuonNoID] = muon.charge();

      if( muon.isGlobalMuon() )     iterL3MuonNoID_isGLB_[_nIterL3MuonNoID] = 1;
      if( muon.isStandAloneMuon() ) iterL3MuonNoID_isSTA_[_nIterL3MuonNoID] = 1;
      if( muon.isTrackerMuon() )    iterL3MuonNoID_isTRK_[_nIterL3MuonNoID] = 1;

      reco::TrackRef innerTrk = muon.innerTrack();
      if( innerTrk.isNonnull() ) {
        iterL3MuonNoID_innerPt_[_nIterL3MuonNoID] = innerTrk->pt();

        tmpTrk trkTmp(innerTrk);
        iterL3NoIDpassed.push_back(trkTmp);

        // const PTrajectoryStateOnDet tmpseed = innerTrk->seedRef()->startingState();
        // tmpTSOD tsod(tmpseed);
        // MuonIterNoIdSeedMap.insert(make_pair(tsod,i));
      }
      else {
        // cout << "IterL3MuonNoID: innerTrk.isNonnull(): this should never happen" << endl;
        tmpTrk trkTmp( -99999. );  // dummy tmpTrk
        iterL3NoIDpassed.push_back(trkTmp);
      }

      _nIterL3MuonNoID++;
    } // -- end of muon iteration

    nIterL3MuonNoID_ = _nIterL3MuonNoID;
  } // -- if getByToken is valid

  //////////////////////////
  // -- Tracks from each algo -- //
  //////////////////////////
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  edm::Handle<TrackingParticleCollection> TPCollection;

  if(iEvent.getByToken(associatorToken, theAssociator) && iEvent.getByToken(trackingParticleToken, TPCollection)){
    fill_trackTemplate(iEvent, t_hltIterL3OIMuonTrack_,          theAssociator, TPCollection, tracker, hltIterL3OIMuonTrackMap,          TThltIterL3OIMuonTrack);
    fill_trackTemplate(iEvent, t_hltIter0IterL3MuonTrack_,       theAssociator, TPCollection, tracker, hltIter0IterL3MuonTrackMap,       TThltIter0IterL3MuonTrack);
    fill_trackTemplate(iEvent, t_hltIter3IterL3MuonTrack_,       theAssociator, TPCollection, tracker, hltIter3IterL3MuonTrackMap,       TThltIter3IterL3MuonTrack);
    fill_trackTemplate(iEvent, t_hltIter0IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, hltIter0IterL3FromL1MuonTrackMap, TThltIter0IterL3FromL1MuonTrack);
    fill_trackTemplate(iEvent, t_hltIter3IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, hltIter3IterL3FromL1MuonTrackMap, TThltIter3IterL3FromL1MuonTrack);

    fill_trackTemplate(iEvent, t_hltIter2IterL3MuonTrack_,       theAssociator, TPCollection, tracker, mvaHltIter2IterL3MuonPixelSeeds_,       hltIter2IterL3MuonTrackMap,      TThltIter2IterL3MuonTrack);
    fill_trackTemplate(iEvent, t_hltIter2IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, mvaHltIter2IterL3FromL1MuonPixelSeeds_, hltIter2IterL3FromL1MuonTrackMap,TThltIter2IterL3FromL1MuonTrack);
  }

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    bool doIso = false;
    if( i==2 ){ //Iter2FromL2Track
      fill_trackTemplate( iEvent, trackCollectionTokens_.at(i), recoToSimCollectionTokens_.at(i), tracker, mvaHltIter2IterL3MuonPixelSeeds_, trkTemplates_.at(i), doIso );
      fill_tpTemplate(    iEvent,                               simToRecoCollectionTokens_.at(i), tracker, mvaHltIter2IterL3MuonPixelSeeds_, tpTemplates_.at(i)         );
    }
    else if( i==4 ){ //Iter2FromL1Track
      fill_trackTemplate( iEvent, trackCollectionTokens_.at(i), recoToSimCollectionTokens_.at(i), tracker, mvaHltIter2IterL3FromL1MuonPixelSeeds_, trkTemplates_.at(i), doIso );
      fill_tpTemplate(    iEvent,                               simToRecoCollectionTokens_.at(i), tracker, mvaHltIter2IterL3FromL1MuonPixelSeeds_, tpTemplates_.at(i)         );
    }
    else{
      fill_trackTemplate( iEvent, trackCollectionTokens_.at(i), recoToSimCollectionTokens_.at(i), trkTemplates_.at(i), doIso );
      fill_tpTemplate(    iEvent,                               simToRecoCollectionTokens_.at(i), tpTemplates_.at(i)         );
    }
  }


}

/*
void MuonHLTNtupler::Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // TrackerHitAssociator associate(iEvent, trackerHitAssociatorConfig_);
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

  //////////////////////////
  // -- hltIterL3OISeedsFromL2Muons -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIterL3OISeedsFromL2Muons;
  if( iEvent.getByToken( t_hltIterL3OISeedsFromL2Muons_, h_hltIterL3OISeedsFromL2Muons) )
  {
    for( auto i=0U; i<h_hltIterL3OISeedsFromL2Muons->size(); ++i )
    {
      const auto& seed(h_hltIterL3OISeedsFromL2Muons->at(i));

      SThltIterL3OISeedsFromL2Muons->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIterL3OISeedsFromL2Muons->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIterL3OIMuonTrackMap.find(seedTsod);
      SThltIterL3OISeedsFromL2Muons->linktmpL3((where2==hltIterL3OIMuonTrackMap.end()) ? -1 : hltIterL3OIMuonTrackMap[seedTsod]);
      SThltIterL3OISeedsFromL2Muons->matchingL3((where2==hltIterL3OIMuonTrackMap.end()) ? -1 : TThltIterL3OIMuonTrack->matchedIDpassedL3(hltIterL3OIMuonTrackMap[seedTsod]));

      // std::cout << "OI RecHit loop start" << std::endl;
      // for (auto rechit = seed.recHits().first; rechit != seed.recHits().second; ++rechit) {
      //   // std::cout << "OI RecHit is valid : " << rechit->isValid() << std::endl;
      //   auto matched = associate.associateHit(*rechit);
      //   // std::cout << "OI Matched RecHit size = " << matched.size() << std::endl;
      // }
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter0IterL3MuonPixelSeedsFromPixelTracks -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter0IterL3MuonPixelSeedsFromPixelTracks;
  if( iEvent.getByToken( t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_, h_hltIter0IterL3MuonPixelSeedsFromPixelTracks) )
  {
    for( auto i=0U; i<h_hltIter0IterL3MuonPixelSeedsFromPixelTracks->size(); ++i )
    {
      const auto& seed(h_hltIter0IterL3MuonPixelSeedsFromPixelTracks->at(i));

      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter0IterL3MuonTrackMap.find(seedTsod);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->linktmpL3((where2==hltIter0IterL3MuonTrackMap.end()) ? -1 : hltIter0IterL3MuonTrackMap[seedTsod]);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->matchingL3((where2==hltIter0IterL3MuonTrackMap.end()) ? -1 : TThltIter0IterL3MuonTrack->matchedIDpassedL3(hltIter0IterL3MuonTrackMap[seedTsod]));

      // std::cout << "IO RecHit loop start" << std::endl;
      // for (auto rechit = seed.recHits().first; rechit != seed.recHits().second; ++rechit) {
      //   // std::cout << "IO RecHit is valid : " << rechit->isValid() << std::endl;
      //   auto matched = associate.associateHit(*rechit);
      //   // std::cout << "IO Matched RecHit size = " << matched.size() << std::endl;
      // }
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter2IterL3MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter2IterL3MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter2IterL3MuonPixelSeeds_, h_hltIter2IterL3MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter2IterL3MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter2IterL3MuonPixelSeeds->at(i));

      SThltIter2IterL3MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter2IterL3MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter2IterL3MuonTrackMap.find(seedTsod);
      SThltIter2IterL3MuonPixelSeeds->linktmpL3((where2==hltIter2IterL3MuonTrackMap.end()) ? -1 : hltIter2IterL3MuonTrackMap[seedTsod]);
      SThltIter2IterL3MuonPixelSeeds->matchingL3((where2==hltIter2IterL3MuonTrackMap.end()) ? -1 : TThltIter2IterL3MuonTrack->matchedIDpassedL3(hltIter2IterL3MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter3IterL3MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter3IterL3MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter3IterL3MuonPixelSeeds_, h_hltIter3IterL3MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter3IterL3MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter3IterL3MuonPixelSeeds->at(i));

      SThltIter3IterL3MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter3IterL3MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter3IterL3MuonTrackMap.find(seedTsod);
      SThltIter3IterL3MuonPixelSeeds->linktmpL3((where2==hltIter3IterL3MuonTrackMap.end()) ? -1 : hltIter3IterL3MuonTrackMap[seedTsod]);
      SThltIter3IterL3MuonPixelSeeds->matchingL3((where2==hltIter3IterL3MuonTrackMap.end()) ? -1 : TThltIter3IterL3MuonTrack->matchedIDpassedL3(hltIter3IterL3MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks;
  if( iEvent.getByToken( t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_, h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks) )
  {
    for( auto i=0U; i<h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->size(); ++i )
    {
      const auto& seed(h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->at(i));

      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter0IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->linktmpL3((where2==hltIter0IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter0IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->matchingL3((where2==hltIter0IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter0IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter0IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter2IterL3FromL1MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter2IterL3FromL1MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter2IterL3FromL1MuonPixelSeeds_, h_hltIter2IterL3FromL1MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter2IterL3FromL1MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter2IterL3FromL1MuonPixelSeeds->at(i));

      SThltIter2IterL3FromL1MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter2IterL3FromL1MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter2IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter2IterL3FromL1MuonPixelSeeds->linktmpL3((where2==hltIter2IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter2IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter2IterL3FromL1MuonPixelSeeds->matchingL3((where2==hltIter2IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter2IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter2IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter3IterL3FromL1MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter3IterL3FromL1MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter3IterL3FromL1MuonPixelSeeds_, h_hltIter3IterL3FromL1MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter3IterL3FromL1MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter3IterL3FromL1MuonPixelSeeds->at(i));

      SThltIter3IterL3FromL1MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter3IterL3FromL1MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter3IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter3IterL3FromL1MuonPixelSeeds->linktmpL3((where2==hltIter3IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter3IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter3IterL3FromL1MuonPixelSeeds->matchingL3((where2==hltIter3IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter3IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter3IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid
}
*/

void MuonHLTNtupler::Fill_TP( const edm::Event &iEvent, tpTemplate* tpTmp )
{
  edm::Handle<TrackingParticleCollection> TPCollection;
  if( iEvent.getByToken(trackingParticleToken, TPCollection) ) {
    for( auto i=0U; i< TPCollection->size(); ++i) {
      if( abs(TPCollection->at(i).pdgId()) == 13 ) {
        tpTmp->fill( TPCollection->at(i) );
      }
    }
  }
}

void MuonHLTNtupler::fill_trackTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
  edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_,
  edm::Handle<TrackingParticleCollection>& TPCollection_,
  edm::ESHandle<TrackerGeometry>& tracker,
  std::map<tmpTSOD,unsigned int>& trkMap,
  trkTemplate* TTtrack
) {
  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( theToken, trkHandle ) )
  {
    auto recSimColl = theAssociator_->associateRecoToSim(trkHandle,TPCollection_);

    for( unsigned int i = 0; i < trkHandle->size(); i++ )
    {
      TTtrack->fill(trkHandle->at(i), bs);

      int linkNo = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3IDpassed.size(); idxL3passed++) {
        if ( iterL3IDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNo = idxL3passed;
      }
      TTtrack->linkIterL3(linkNo);

      int linkNoNoId = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3NoIDpassed.size(); idxL3passed++) {
        if ( iterL3NoIDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNoNoId = idxL3passed;
      }
      TTtrack->linkIterL3NoId(linkNoNoId);

      const PTrajectoryStateOnDet tmpseed = trkHandle->at(i).seedRef()->startingState();
      tmpTSOD tsod(tmpseed);
      trkMap.insert(make_pair(tsod,i));

      auto track = trkHandle->refAt(i);
      auto TPfound = recSimColl.find(track);
      if (TPfound != recSimColl.end()) {
        const auto& TPmatch = TPfound->val;
        TTtrack->fillBestTP(TPmatch[0].first);
        TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
        TTtrack->fillmatchedTPsize(TPmatch.size());
      } else {  // sync vector size
        TTtrack->fillDummyTP();
        TTtrack->fillBestTPsharedFrac(-99999.);
        TTtrack->fillmatchedTPsize(0);
      }
    }
  }
}

void MuonHLTNtupler::fill_trackTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
  edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_,
  edm::Handle<TrackingParticleCollection>& TPCollection_,
  edm::ESHandle<TrackerGeometry>& tracker,
  pairSeedMvaEstimator pairMvaEstimator,
  std::map<tmpTSOD,unsigned int>& trkMap,
  trkTemplate* TTtrack
) {

  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);
  const l1t::MuonBxCollection l1Muons = *(h_L1Muon.product());

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );
  const reco::RecoChargedCandidateCollection l2Muons = *(h_L2Muon.product());

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( theToken, trkHandle ) )
  {
    auto recSimColl = theAssociator_->associateRecoToSim(trkHandle,TPCollection_);

    for( unsigned int i = 0; i < trkHandle->size(); i++ )
    {
      TTtrack->fill(trkHandle->at(i), bs);

      int linkNo = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3IDpassed.size(); idxL3passed++) {
        if ( iterL3IDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNo = idxL3passed;
      }
      TTtrack->linkIterL3(linkNo);

      int linkNoNoId = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3NoIDpassed.size(); idxL3passed++) {
        if ( iterL3NoIDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNoNoId = idxL3passed;
      }
      TTtrack->linkIterL3NoId(linkNoNoId);

      const PTrajectoryStateOnDet tmpseed = trkHandle->at(i).seedRef()->startingState();
      tmpTSOD tsod(tmpseed);
      trkMap.insert(make_pair(tsod,i));

      auto track = trkHandle->refAt(i);
      auto TPfound = recSimColl.find(track);
      if (TPfound != recSimColl.end()) {
        const auto& TPmatch = TPfound->val;
        TTtrack->fillBestTP(TPmatch[0].first);
        TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
        TTtrack->fillmatchedTPsize(TPmatch.size());
      } else {  // sync vector size
        TTtrack->fillDummyTP();
        TTtrack->fillBestTPsharedFrac(-99999.);
        TTtrack->fillmatchedTPsize(0);
      }

      if( doMVA && hasL1 && hasL2 ) {
        const TrajectorySeed seed = *(trkHandle->at(i).seedRef());
        GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
        //GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

        vector<double> mva = getSeedMva(
          pairMvaEstimator,
          seed,
          global_p,
          l1Muons,
          l2Muons
        );
        //TTtrack->fillMva( mva[0], mva[1], mva[2], mva[3] );
        TTtrack->fillMva( mva[0], -99999., -99999., -99999. );
      }
      else {
        // cout << "fill_trackTemplate: !(hasL1 && hasL2 && hasL1TkMu)" << endl;
        TTtrack->fillMva( -99999., -99999., -99999., -99999. );
      }
    }
  }
}

void MuonHLTNtupler::fill_trackTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& trkToken,
  edm::EDGetTokenT<reco::RecoToSimCollection>& assoToken,
  trkTemplate* TTtrack,
  bool doIso = false
) {

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( trkToken, trkHandle ) ) {

    edm::Handle<reco::RecoToSimCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto recSimColl = *assoHandle.product();

      for( unsigned int i = 0; i < trkHandle->size(); i++ ) {
        TTtrack->fill(trkHandle->at(i), bs);

        auto track = trkHandle->refAt(i);
        auto TPfound = recSimColl.find(track);
        if (TPfound != recSimColl.end()) {
          const auto& TPmatch = TPfound->val;
          TTtrack->fillBestTP(TPmatch[0].first);
          TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
          TTtrack->fillmatchedTPsize(TPmatch.size());
        } else {  // to sync vector size
          TTtrack->fillDummyTP();
          TTtrack->fillBestTPsharedFrac(-99999.);
          TTtrack->fillmatchedTPsize(0);
        }

        // -- fill dummy
        TTtrack->linkIterL3(-1);
        TTtrack->linkIterL3NoId(-1);
        TTtrack->fillMva( -99999., -99999., -99999., -99999. );
      }
    }
  }
}

void MuonHLTNtupler::fill_trackTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& trkToken,
  edm::EDGetTokenT<reco::RecoToSimCollection>& assoToken,
  edm::ESHandle<TrackerGeometry>& tracker,
  pairSeedMvaEstimator pairMvaEstimator,
  trkTemplate* TTtrack,
  bool doIso = false
) {

  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);
  const l1t::MuonBxCollection l1Muons = *(h_L1Muon.product());

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );
  const reco::RecoChargedCandidateCollection l2Muons = *(h_L2Muon.product());

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( trkToken, trkHandle ) ) {

    edm::Handle<reco::RecoToSimCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto recSimColl = *assoHandle.product();

      for( unsigned int i = 0; i < trkHandle->size(); i++ ) {
        TTtrack->fill(trkHandle->at(i), bs);

        // -- fill dummy index
        TTtrack->linkIterL3(-1);
        TTtrack->linkIterL3NoId(-1);

        auto track = trkHandle->refAt(i);
        auto TPfound = recSimColl.find(track);
        if (TPfound != recSimColl.end()) {
          const auto& TPmatch = TPfound->val;
          TTtrack->fillBestTP(TPmatch[0].first);
          TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
          TTtrack->fillmatchedTPsize(TPmatch.size());
        } else {  // to sync vector size
          TTtrack->fillDummyTP();
          TTtrack->fillBestTPsharedFrac(-99999.);
          TTtrack->fillmatchedTPsize(0);
        }

        if( doMVA && hasL1 && hasL2 ) {
          const TrajectorySeed seed = *(trkHandle->at(i).seedRef());
          GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
          //GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());
  
          vector<double> mva = getSeedMva(
            pairMvaEstimator,
            seed,
            global_p,
            l1Muons,
            l2Muons
          );
          //TTtrack->fillMva( mva[0], mva[1], mva[2], mva[3] );
          TTtrack->fillMva( mva[0], -99999., -99999., -99999. );
        }
        else {
          // cout << "fill_trackTemplate: !(hasL1 && hasL2 && hasL1TkMu)" << endl;
          TTtrack->fillMva( -99999., -99999., -99999., -99999. );
        }
      }
    }
  }
}

void MuonHLTNtupler::fill_tpTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<reco::SimToRecoCollection>& assoToken,
  tpTemplate* TTtp
) {

  edm::Handle<TrackingParticleCollection> TPCollection;
  if( iEvent.getByToken(trackingParticleToken, TPCollection) ) {

    edm::Handle<reco::SimToRecoCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto simRecColl = *assoHandle.product();

      for( unsigned int i = 0; i < TPCollection->size(); i++ ) {

        auto tp = TPCollection->at(i);

        if( !(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0) )
          continue;

        if( abs( tp.pdgId() ) != 13 )
          continue;

        bool isStable = true;
        for (TrackingParticle::genp_iterator j = tp.genParticle_begin(); j != tp.genParticle_end(); ++j) {
          if (j->get() == nullptr || j->get()->status() != 1) {
            isStable = false;
          }
        }
        if( tp.status() == -99 && (std::abs(tp.pdgId()) != 11 && std::abs(tp.pdgId()) != 13 && std::abs(tp.pdgId()) != 211 &&
                                   std::abs(tp.pdgId()) != 321 && std::abs(tp.pdgId()) != 2212 && std::abs(tp.pdgId()) != 3112 &&
                                   std::abs(tp.pdgId()) != 3222 && std::abs(tp.pdgId()) != 3312 && std::abs(tp.pdgId()) != 3334)
        ) {
          isStable = false;
        }

        if( !isStable )
          continue;

        TTtp->fill( tp );

        TrackingParticleRef tpref(TPCollection, i);
        auto TrkFound = simRecColl.find( tpref );
        if( TrkFound != simRecColl.end() ) {
          const auto& trkMatch = TrkFound->val;
          TTtp->fill_matchedTrk(
            trkMatch[0].first->pt(),
            trkMatch[0].first->eta(),
            trkMatch[0].first->phi(),
            trkMatch[0].first->charge(),
            trkMatch[0].first->px(),
            trkMatch[0].first->py(),
            trkMatch[0].first->pz(),
            trkMatch[0].first->vx(),
            trkMatch[0].first->vy(),
            trkMatch[0].first->vz(),
            trkMatch[0].first->dxy(bs->position()),
            trkMatch[0].first->dxyError(*bs),
            trkMatch[0].first->dz(bs->position()),
            trkMatch[0].first->dzError(),
            trkMatch[0].first->normalizedChi2(),
            trkMatch[0].second,
            trkMatch[0].first->numberOfValidHits(),
            -99999.
          );
        }
        else {
          TTtp->fill_matchedTrk(
            -99999.,
            -99999.,
            -99999.,
            -99999,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999,
            -99999.
          );
        }
      }
    }
  }
}

void MuonHLTNtupler::fill_tpTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<reco::SimToRecoCollection>& assoToken,
  edm::ESHandle<TrackerGeometry>& tracker,
  pairSeedMvaEstimator pairMvaEstimator,
  tpTemplate* TTtp
) {

  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);
  const l1t::MuonBxCollection l1Muons = *(h_L1Muon.product());

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );
  const reco::RecoChargedCandidateCollection l2Muons = *(h_L2Muon.product());

  edm::Handle<TrackingParticleCollection> TPCollection;
  if( iEvent.getByToken(trackingParticleToken, TPCollection) ) {

    edm::Handle<reco::SimToRecoCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto simRecColl = *assoHandle.product();

      for( unsigned int i = 0; i < TPCollection->size(); i++ ) {

        auto tp = TPCollection->at(i);

        if( !(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0) )
          continue;

        if( abs( tp.pdgId() ) != 13 )
          continue;

        bool isStable = true;
        for (TrackingParticle::genp_iterator j = tp.genParticle_begin(); j != tp.genParticle_end(); ++j) {
          if (j->get() == nullptr || j->get()->status() != 1) {
            isStable = false;
          }
        }
        if( tp.status() == -99 && (std::abs(tp.pdgId()) != 11 && std::abs(tp.pdgId()) != 13 && std::abs(tp.pdgId()) != 211 &&
                                   std::abs(tp.pdgId()) != 321 && std::abs(tp.pdgId()) != 2212 && std::abs(tp.pdgId()) != 3112 &&
                                   std::abs(tp.pdgId()) != 3222 && std::abs(tp.pdgId()) != 3312 && std::abs(tp.pdgId()) != 3334)
        ) {
          isStable = false;
        }

        if( !isStable )
          continue;

        TTtp->fill( tp );

        TrackingParticleRef tpref(TPCollection, i);
        auto TrkFound = simRecColl.find( tpref );
        if( TrkFound != simRecColl.end() ) {
          const auto& trkMatch = TrkFound->val;

          double this_mva = -99999.;
          if( doMVA && hasL1 && hasL2 ) {
            const TrajectorySeed seed = *(trkMatch[0].first->seedRef());
            GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
            //GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());
    
            vector<double> mva = getSeedMva(
              pairMvaEstimator,
              seed,
              global_p,
              l1Muons,
              l2Muons
            );
            this_mva = mva[0];
          }

          TTtp->fill_matchedTrk(
            trkMatch[0].first->pt(),
            trkMatch[0].first->eta(),
            trkMatch[0].first->phi(),
            trkMatch[0].first->charge(),
            trkMatch[0].first->px(),
            trkMatch[0].first->py(),
            trkMatch[0].first->pz(),
            trkMatch[0].first->vx(),
            trkMatch[0].first->vy(),
            trkMatch[0].first->vz(),
            trkMatch[0].first->dxy(bs->position()),
            trkMatch[0].first->dxyError(*bs),
            trkMatch[0].first->dz(bs->position()),
            trkMatch[0].first->dzError(),
            trkMatch[0].first->normalizedChi2(),
            trkMatch[0].second,
            trkMatch[0].first->numberOfValidHits(),
            this_mva
          );
        }
        else {
          TTtp->fill_matchedTrk(
            -99999.,
            -99999.,
            -99999.,
            -99999,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999.,
            -99999,
            -99999.
          );
        }
      }
    }
  }
}

// -- reference: https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/src/MuonSelectors.cc#L910-L938
bool MuonHLTNtupler::isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx){
  if(!muon.isGlobalMuon()) return false;

  bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                     muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );

  bool muMatchedSt = muon.numberOfMatchedStations()>1;
  if(!muMatchedSt) {
    if( muon.isTrackerMuon() && muon.numberOfMatchedStations()==1 ) {
      if( muon.expectedNnumberOfMatchedStations()<2 ||
          !(muon.stationMask()==1 || muon.stationMask()==16) ||
          muon.numberOfMatchedRPCLayers()>2
        )
        muMatchedSt = true;
    }
  }

  bool muID = muValHits && muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

void MuonHLTNtupler::endJob() {
  //for( int i=0; i<4; ++i ) {
  for( int i=0; i<1; ++i ) {
    delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).first;
    delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).second;
    delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).first;
    delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).second;
  }

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    delete trkTemplates_.at(i);
    delete tpTemplates_.at(i);
  }
}

void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
