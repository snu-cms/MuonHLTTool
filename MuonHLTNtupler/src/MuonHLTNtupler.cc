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
// SaveAllTracks(iConfig.getParameter<bool>("SaveAllTracks")),
// SaveStubs(iConfig.getParameter<bool>("SaveStubs")),
// ttTrackToken_        ( consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"     )) ),
// ttTrackMCTruthToken_ ( consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag"))),
// ttStubToken_         ( consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter<edm::InputTag>("L1StubInputTag"))),
// TkMuonToken_         ( consumes<l1t::TkMuonCollection>          (iConfig.getParameter<edm::InputTag>("TkMuonToken"))),
// l1TkPrimaryVertexToken_ (consumes< l1t::TkPrimaryVertexCollection> (iConfig.getParameter<edm::InputTag>("l1TkPrimaryVertex"))),

trackerHitAssociatorConfig_(iConfig, consumesCollector()),
associatorToken(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("associator"))),
trackingParticleToken(consumes<TrackingParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingParticle"))),
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
t_genParticle_       ( consumes< reco::GenParticleCollection >            (iConfig.getUntrackedParameter<edm::InputTag>("genParticle"       )) )
{
  /*
  mvaFileHltIterL3OISeedsFromL2Muons_B_0_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_B_0");
  mvaFileHltIterL3OISeedsFromL2Muons_B_1_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_B_1");
  mvaFileHltIterL3OISeedsFromL2Muons_B_2_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_B_2");
  mvaFileHltIterL3OISeedsFromL2Muons_B_3_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_B_3");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3");
  mvaFileHltIter2IterL3MuonPixelSeeds_B_0_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_B_0");
  mvaFileHltIter2IterL3MuonPixelSeeds_B_1_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_B_1");
  mvaFileHltIter2IterL3MuonPixelSeeds_B_2_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_B_2");
  mvaFileHltIter2IterL3MuonPixelSeeds_B_3_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_B_3");
  mvaFileHltIter3IterL3MuonPixelSeeds_B_0_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_B_0");
  mvaFileHltIter3IterL3MuonPixelSeeds_B_1_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_B_1");
  mvaFileHltIter3IterL3MuonPixelSeeds_B_2_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_B_2");
  mvaFileHltIter3IterL3MuonPixelSeeds_B_3_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_B_3");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3");
  mvaFileHltIterL3OISeedsFromL2Muons_E_0_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_E_0");
  mvaFileHltIterL3OISeedsFromL2Muons_E_1_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_E_1");
  mvaFileHltIterL3OISeedsFromL2Muons_E_2_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_E_2");
  mvaFileHltIterL3OISeedsFromL2Muons_E_3_                       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIterL3OISeedsFromL2Muons_E_3");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2");
  mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3_       = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3");
  mvaFileHltIter2IterL3MuonPixelSeeds_E_0_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_E_0");
  mvaFileHltIter2IterL3MuonPixelSeeds_E_1_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_E_1");
  mvaFileHltIter2IterL3MuonPixelSeeds_E_2_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_E_2");
  mvaFileHltIter2IterL3MuonPixelSeeds_E_3_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3MuonPixelSeeds_E_3");
  mvaFileHltIter3IterL3MuonPixelSeeds_E_0_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_E_0");
  mvaFileHltIter3IterL3MuonPixelSeeds_E_1_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_E_1");
  mvaFileHltIter3IterL3MuonPixelSeeds_E_2_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_E_2");
  mvaFileHltIter3IterL3MuonPixelSeeds_E_3_                      = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3MuonPixelSeeds_E_3");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2");
  mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3_ = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2");
  mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3_                = iConfig.getParameter<edm::FileInPath>("mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3");

  mvaScaleMeanHltIterL3OISeedsFromL2Muons_B_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIterL3OISeedsFromL2Muons_B");
  mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_ =       iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B");
  mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B");
  mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B");
  mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_ = iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B");
  mvaScaleStdHltIterL3OISeedsFromL2Muons_B_ =                        iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIterL3OISeedsFromL2Muons_B");
  mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_ =        iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B");
  mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3MuonPixelSeeds_B");
  mvaScaleStdHltIter3IterL3MuonPixelSeeds_B_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter3IterL3MuonPixelSeeds_B");
  mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_ =  iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B");
  mvaScaleMeanHltIterL3OISeedsFromL2Muons_E_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIterL3OISeedsFromL2Muons_E");
  mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_ =       iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E");
  mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E");
  mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E_ =                      iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E");
  mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_ = iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E");
  mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E_ =                iConfig.getParameter<std::vector<double>>("mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E");
  mvaScaleStdHltIterL3OISeedsFromL2Muons_E_ =                        iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIterL3OISeedsFromL2Muons_E");
  mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_ =        iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E");
  mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3MuonPixelSeeds_E");
  mvaScaleStdHltIter3IterL3MuonPixelSeeds_E_ =                       iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter3IterL3MuonPixelSeeds_E");
  mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_ =  iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E");
  mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E_ =                 iConfig.getParameter<std::vector<double>>("mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E");
  */
}

void MuonHLTNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init();

  // -- basic info.
  isRealData_ = iEvent.isRealData();

  runNum_       = iEvent.id().run();
  lumiBlockNum_ = iEvent.id().luminosityBlock();
  eventNum_     = iEvent.id().event();

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
  // Fill_HLT(iEvent, 1); // -- rerun objects
  Fill_HLTMuon(iEvent);
  Fill_L1Muon(iEvent);
  Fill_IterL3(iEvent, iSetup);
  if( doSeed )  Fill_Seed(iEvent, iSetup);
  if( !isRealData_ ) {
    Fill_GenParticle(iEvent);
    Fill_TP(iEvent, TrkParticle);
  }

  ntuple_->Fill();
}

void MuonHLTNtupler::beginJob()
{
  /*
  mvaHltIterL3OISeedsFromL2Muons_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_B_0_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_B_, mvaScaleStdHltIterL3OISeedsFromL2Muons_B_),
               new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_E_0_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_E_, mvaScaleStdHltIterL3OISeedsFromL2Muons_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_B_1_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_B_, mvaScaleStdHltIterL3OISeedsFromL2Muons_B_),
               new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_E_1_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_E_, mvaScaleStdHltIterL3OISeedsFromL2Muons_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_B_2_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_B_, mvaScaleStdHltIterL3OISeedsFromL2Muons_B_),
               new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_E_2_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_E_, mvaScaleStdHltIterL3OISeedsFromL2Muons_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_B_3_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_B_, mvaScaleStdHltIterL3OISeedsFromL2Muons_B_),
               new SeedMvaEstimator(mvaFileHltIterL3OISeedsFromL2Muons_E_3_, mvaScaleMeanHltIterL3OISeedsFromL2Muons_E_, mvaScaleStdHltIterL3OISeedsFromL2Muons_E_) )
  };
  mvaHltIter0IterL3MuonPixelSeedsFromPixelTracks_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3_, mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_) )
  };
  mvaHltIter2IterL3MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_B_0_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_E_0_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_B_1_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_E_1_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_B_2_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_E_2_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_B_3_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3MuonPixelSeeds_E_3_, mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3MuonPixelSeeds_E_) )
  };
  mvaHltIter3IterL3MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_B_0_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_E_0_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_B_1_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_E_1_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_B_2_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_E_2_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_B_3_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3MuonPixelSeeds_E_3_, mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3MuonPixelSeeds_E_) )
  };
  mvaHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_),
               new SeedMvaEstimator(mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3_, mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_, mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_) )
  };
  mvaHltIter2IterL3FromL1MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_) )
  };
  mvaHltIter3IterL3FromL1MuonPixelSeeds_ = {
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E_) ),
    make_pair( new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B_),
               new SeedMvaEstimator(mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3_, mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E_) )
  };
  */

  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");

  // initilize
  /*
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_bendchi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_lhits = new std::vector<int>;
  m_trk_dhits = new std::vector<int>;
  m_trk_seed    = new std::vector<int>;
  m_trk_phiSector    = new std::vector<unsigned int>;
  m_trk_genuine       = new std::vector<int>;
  m_trk_loose         = new std::vector<int>;
  m_trk_unknown       = new std::vector<int>;
  m_trk_combinatoric  = new std::vector<int>;
  m_trk_fake          = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt    = new std::vector<float>;
  m_trk_matchtp_eta   = new std::vector<float>;
  m_trk_matchtp_phi   = new std::vector<float>;
  m_trk_matchtp_z0    = new std::vector<float>;
  m_trk_matchtp_dxy   = new std::vector<float>;

  m_stub_x = new std::vector<float>;
  m_stub_y = new std::vector<float>;
  m_stub_z = new std::vector<float>;
  m_stub_isBarrel = new std::vector<int>;
  m_stub_layer    = new std::vector<int>;
  */

  Make_Branch();
}

void MuonHLTNtupler::Init()
{
  // if (SaveAllTracks){
  // m_trk_pt->clear();
  // m_trk_eta->clear();
  // m_trk_phi->clear();
  // m_trk_d0->clear();
  // m_trk_z0->clear();
  // m_trk_rInv.clear();
  // m_trk_tanL.clear();
  // m_trk_MVA1.clear();
  // m_trk_MVA2.clear();
  // m_trk_MVA3.clear();
  // m_trk_chi2->clear();
  // m_trk_bendchi2->clear();
  // m_trk_nstub->clear();
  // m_trk_lhits->clear();
  // m_trk_dhits->clear();
  // m_trk_seed->clear();
  // m_trk_phiSector->clear();
  // m_trk_genuine->clear();
  // m_trk_loose->clear();
  // m_trk_unknown->clear();
  // m_trk_combinatoric->clear();
  // m_trk_fake->clear();
  // m_trk_matchtp_pdgid->clear();
  // m_trk_matchtp_pt->clear();
  // m_trk_matchtp_eta->clear();
  // m_trk_matchtp_phi->clear();
  // m_trk_matchtp_z0->clear();
  // m_trk_matchtp_dxy->clear();

  // m_stub_x->clear();
  // m_stub_y->clear();
  // m_stub_z->clear();
  // m_stub_isBarrel->clear();
  // m_stub_layer->clear();

  // // l1TkMuon
  // mL1TkMu_pt.clear();
  // mL1TkMu_eta.clear();
  // mL1TkMu_phi.clear();
  // mL1TkMu_trkIsol.clear();
  // mL1TkMu_trkzVtx.clear();
  // mL1TkMu_dR.clear();

  // mL1TkMu_nTracksMatched.clear();
  // mL1TkMu_trackCurvature.clear();

  // mL1TkMu_quality.clear();
  // mL1TkMu_pattern.clear();
  // mL1TkMu_muonDetector.clear();

  // mL1TkMu_TTTpointer.clear();

  // mL1TkMu_muRefHwPt.clear();
  // mL1TkMu_muRefHwDXY.clear();
  // mL1TkMu_muRefHwEta.clear();
  // mL1TkMu_muRefHwPhi.clear();
  // mL1TkMu_muRefHwSign.clear();
  // mL1TkMu_muRefHwSignValid.clear();
  // mL1TkMu_muRefHwQual.clear();

  // mTTTrackMap.clear();
  // }


  isRealData_ = false;

  runNum_       = -999;
  lumiBlockNum_ = -999;
  eventNum_     = 0;

  bunchID_ = -999;

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
  }

  // -- original trigger objects- - //
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
  // MuonIterNoIdSeedMap.clear();
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

  SThltIterL3OISeedsFromL2Muons->clear();
  SThltIter0IterL3MuonPixelSeedsFromPixelTracks->clear();
  SThltIter2IterL3MuonPixelSeeds->clear();
  SThltIter3IterL3MuonPixelSeeds->clear();
  SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->clear();
  SThltIter2IterL3FromL1MuonPixelSeeds->clear();
  SThltIter3IterL3FromL1MuonPixelSeeds->clear();

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
}

void MuonHLTNtupler::Make_Branch()
{
  /*
  if (SaveAllTracks){
  ntuple_->Branch("trk_pt",    &m_trk_pt);
  ntuple_->Branch("trk_eta",   &m_trk_eta);
  ntuple_->Branch("trk_phi",   &m_trk_phi);
  ntuple_->Branch("trk_d0",    &m_trk_d0);
  ntuple_->Branch("trk_z0",    &m_trk_z0);
  ntuple_->Branch("trk_rInv", &m_trk_rInv);
  ntuple_->Branch("trk_tanL", &m_trk_tanL);
  ntuple_->Branch("trk_MVA1", &m_trk_MVA1);
  ntuple_->Branch("trk_MVA2", &m_trk_MVA2);
  ntuple_->Branch("trk_MVA3", &m_trk_MVA3);

  ntuple_->Branch("trk_chi2",  &m_trk_chi2);
  ntuple_->Branch("trk_bendchi2",  &m_trk_bendchi2);
  ntuple_->Branch("trk_nstub", &m_trk_nstub);
  ntuple_->Branch("trk_lhits", &m_trk_lhits);
  ntuple_->Branch("trk_dhits", &m_trk_dhits);
  ntuple_->Branch("trk_seed",    &m_trk_seed);
  ntuple_->Branch("trk_phiSector", &m_trk_phiSector);
  ntuple_->Branch("trk_genuine",      &m_trk_genuine);
  ntuple_->Branch("trk_loose",        &m_trk_loose);
  ntuple_->Branch("trk_unknown",      &m_trk_unknown);
  ntuple_->Branch("trk_combinatoric", &m_trk_combinatoric);
  ntuple_->Branch("trk_fake",         &m_trk_fake);
  ntuple_->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
  ntuple_->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
  ntuple_->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
  ntuple_->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
  ntuple_->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
  ntuple_->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);

  ntuple_->Branch("stub_x", &m_stub_x);
  ntuple_->Branch("stub_y", &m_stub_y);
  ntuple_->Branch("stub_z", &m_stub_z);
  ntuple_->Branch("stub_isBarrel",   &m_stub_isBarrel);
  ntuple_->Branch("stub_layer",      &m_stub_layer);

  ntuple_->Branch("L1TkMu_pt", &mL1TkMu_pt);
  ntuple_->Branch("L1TkMu_eta", &mL1TkMu_eta);
  ntuple_->Branch("L1TkMu_phi", &mL1TkMu_phi);

  ntuple_->Branch("L1TkMu_trkIsol", &mL1TkMu_trkIsol);
  ntuple_->Branch("L1TkMu_trkzVtx", &mL1TkMu_trkzVtx);
  ntuple_->Branch("L1TkMu_dR", &mL1TkMu_dR);

  ntuple_->Branch("L1TkMu_nTracksMatched", &mL1TkMu_nTracksMatched);
  ntuple_->Branch("L1TkMu_trackCurvature", &mL1TkMu_trackCurvature);

  ntuple_->Branch("L1TkMu_quality", &mL1TkMu_quality);
  ntuple_->Branch("L1TkMu_pattern", &mL1TkMu_pattern);
  ntuple_->Branch("L1TkMu_muonDetector", &mL1TkMu_muonDetector);

  ntuple_->Branch("L1TkMu_TTTpointer", &mL1TkMu_TTTpointer);

  ntuple_->Branch("L1TkMu_muRefHwPt", &mL1TkMu_muRefHwPt);
  ntuple_->Branch("L1TkMu_muRefHwDXY", &mL1TkMu_muRefHwDXY);
  ntuple_->Branch("L1TkMu_muRefHwEta", &mL1TkMu_muRefHwEta);
  ntuple_->Branch("L1TkMu_muRefHwPhi", &mL1TkMu_muRefHwPhi);
  ntuple_->Branch("L1TkMu_muRefHwSign", &mL1TkMu_muRefHwSign);
  ntuple_->Branch("L1TkMu_muRefHwSignValid", &mL1TkMu_muRefHwSignValid);
  ntuple_->Branch("L1TkMu_muRefHwQual", &mL1TkMu_muRefHwQual);
  }
  */

  ntuple_->Branch("isRealData", &isRealData_, "isRealData/O"); // -- O: boolean -- //
  ntuple_->Branch("runNum",&runNum_,"runNum/I");
  ntuple_->Branch("lumiBlockNum",&lumiBlockNum_,"lumiBlockNum/I");
  ntuple_->Branch("eventNum",&eventNum_,"eventNum/l"); // -- unsigned long long -- //
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

  SThltIterL3OISeedsFromL2Muons->setBranch(ntuple_,"hltIterL3OISeedsFromL2Muons");
  SThltIter0IterL3MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3MuonPixelSeedsFromPixelTracks");
  SThltIter2IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3MuonPixelSeeds");
  SThltIter3IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3MuonPixelSeeds");
  SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks");
  SThltIter2IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3FromL1MuonPixelSeeds");
  SThltIter3IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3FromL1MuonPixelSeeds");

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
}

/*
void MuonHLTNtupler::Fill_L1Track(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  if (SaveStubs) iEvent.getByToken(ttStubToken_, TTStubHandle);

  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);
  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  if (SaveAllTracks){
  if (DebugMode) {
    cout << endl << "Loop over L1 tracks!" << endl;
    // cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
  }

  int this_l1track = 0;
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
  for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);

    float tmp_trk_pt   = iterL1Track->momentum().perp();
    float tmp_trk_eta  = iterL1Track->momentum().eta();
    float tmp_trk_phi  = iterL1Track->momentum().phi();
    float tmp_trk_z0   = iterL1Track->z0(); //cm

    float tmp_trk_rInv = static_cast<float>(iterL1Track->rInv());
    float tmp_trk_tanL = static_cast<float>(iterL1Track->tanL());
    float tmp_trk_MVA1 = static_cast<float>(iterL1Track->trkMVA1());
    float tmp_trk_MVA2 = static_cast<float>(iterL1Track->trkMVA2());
    float tmp_trk_MVA3 = static_cast<float>(iterL1Track->trkMVA3());

    float tmp_trk_d0 = static_cast<float>(iterL1Track->d0());

    // if (L1Tk_nPar == 5) {
    //   float tmp_trk_x0   = iterL1Track->POCA().x();
    //   float tmp_trk_y0   = iterL1Track->POCA().y();
    //   tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
    //   }

    float tmp_trk_chi2 = iterL1Track->chi2();
    float tmp_trk_bendchi2 = iterL1Track->stubPtConsistency();

    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub  = (int) stubRefs.size();

    int tmp_trk_seed = 0;
    // if (SaveTracklet) tmp_trk_seed = (int) iterL1Track->trackSeedType();

    unsigned int tmp_trk_phiSector = iterL1Track->phiSector();

    // ----------------------------------------------------------------------------------------------
    // loop over stubs on tracks

    //float tmp_trk_bend_chi2 = 0;
    int tmp_trk_dhits=0;
    int tmp_trk_lhits=0;

    // loop over stubs
    for (int is=0; is<tmp_trk_nstub; is++) {

      //detID of stub
      DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->clusterRef(0))->getDetId() )->geographicalId();

      MeasurementPoint coords = stubRefs.at(is)->clusterRef(0)->findAverageLocalCoordinatesCentered();
      const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
      Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );

      double x=posStub.x();
      double y=posStub.y();
      double z=posStub.z();

      int isBarrel = 0;
      int layer=-999999;
      if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
        isBarrel = 1;
        layer  = static_cast<int>(tTopo->layer(detIdStub));
        if (DebugMode) cout << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
        tmp_trk_lhits+=pow(10,layer-1);
      }
      else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
        isBarrel = 0;
        layer  = static_cast<int>(tTopo->layer(detIdStub));
        if (DebugMode) cout << "   stub in disk " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
        tmp_trk_dhits+=pow(10,layer-1);
      }

      m_stub_x->push_back(x);
      m_stub_y->push_back(y);
      m_stub_z->push_back(z);
      m_stub_isBarrel->push_back(isBarrel);
      m_stub_layer->push_back(layer);
    }//end loop over stubs

    int tmp_trk_genuine = 0;
    int tmp_trk_loose = 0;
    int tmp_trk_unknown = 0;
    int tmp_trk_combinatoric = 0;
    if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
    if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
    if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
    if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;

    if (DebugMode) {
    cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
       << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
    if (tmp_trk_genuine) cout << " (is genuine)" << endl;
    if (tmp_trk_unknown) cout << " (is unknown)" << endl;
    if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl;
    }

  m_trk_pt ->push_back(tmp_trk_pt);
  m_trk_eta->push_back(tmp_trk_eta);
  m_trk_phi->push_back(tmp_trk_phi);
  m_trk_z0 ->push_back(tmp_trk_z0);
  m_trk_d0->push_back(tmp_trk_d0);
  m_trk_rInv.push_back(tmp_trk_rInv);
  m_trk_tanL.push_back(tmp_trk_tanL);
  m_trk_MVA1.push_back(tmp_trk_MVA1);
  m_trk_MVA2.push_back(tmp_trk_MVA2);
  m_trk_MVA3.push_back(tmp_trk_MVA3);
  m_trk_chi2 ->push_back(tmp_trk_chi2);
  m_trk_bendchi2 ->push_back(tmp_trk_bendchi2);
  m_trk_nstub->push_back(tmp_trk_nstub);
  // m_trk_dhits->push_back(tmp_trk_dhits);
  // m_trk_lhits->push_back(tmp_trk_lhits);
  m_trk_seed->push_back(tmp_trk_seed);
  m_trk_phiSector->push_back(tmp_trk_phiSector);
  m_trk_genuine->push_back(tmp_trk_genuine);
  m_trk_loose->push_back(tmp_trk_loose);
  m_trk_unknown->push_back(tmp_trk_unknown);
  m_trk_combinatoric->push_back(tmp_trk_combinatoric);

  MuonHLTobjCorrelator::L1TTTrack theTrack( *iterL1Track );

  mTTTrackMap.insert( make_pair(theTrack, static_cast<unsigned int>(this_l1track)) );

  // ----------------------------------------------------------------------------------------------
  // for studying the fake rate
  // ----------------------------------------------------------------------------------------------

  edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

  int myFake = 0;

  int myTP_pdgid = -999;
  float myTP_pt = -999;
  float myTP_eta = -999;
  float myTP_phi = -999;
  float myTP_z0 = -999;
  float myTP_dxy = -999;

  if (my_tp.isNull()) myFake = 0;
  else {
  int tmp_eventid = my_tp->eventId().event();

  if (tmp_eventid > 0) myFake = 2;
  else myFake = 1;

  myTP_pdgid = my_tp->pdgId();
  myTP_pt = my_tp->p4().pt();
  myTP_eta = my_tp->p4().eta();
  myTP_phi = my_tp->p4().phi();
  myTP_z0 = my_tp->vertex().z();

  float myTP_x0 = my_tp->vertex().x();
  float myTP_y0 = my_tp->vertex().y();
  myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

  if (DebugMode) {
    cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta()
         << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
         << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
        }
      }

  m_trk_fake->push_back(myFake);

  m_trk_matchtp_pdgid->push_back(myTP_pdgid);
  m_trk_matchtp_pt->push_back(myTP_pt);
  m_trk_matchtp_eta->push_back(myTP_eta);
  m_trk_matchtp_phi->push_back(myTP_phi);
  m_trk_matchtp_z0->push_back(myTP_z0);
  m_trk_matchtp_dxy->push_back(myTP_dxy);

  this_l1track++;
    }//l1 track loop
  }//save all tracks

  edm::Handle<l1t::TkMuonCollection> TkMuon;
  iEvent.getByToken(TkMuonToken_,TkMuon);
  for(auto Tkmu=TkMuon->begin(); Tkmu!=TkMuon->end(); ++Tkmu)
  {
    // https://github.com/cms-sw/cmssw/blob/09b17fcfb3900782ab78ad6e0c76e1957c94ff71/DataFormats/L1TCorrelator/interface/TkMuon.h
    // https://github.com/cms-sw/cmssw/blob/2ee90040994dde0e6a4fe686de614b09916bd80a/DataFormats/L1TrackTrigger/interface/TTTrack.h
    // cout<< "phi"<<Tkmu->trkPtr()->phi()<<endl;
    // cout<< "eta"<<Tkmu->trkPtr()->eta()<<endl;
    // cout<< "z0"<<Tkmu->trkPtr()->z0()<<endl;

    // l1TkMuon
    mL1TkMu_pt.push_back( Tkmu->pt() );
    mL1TkMu_eta.push_back( Tkmu->eta() );
    mL1TkMu_phi.push_back( Tkmu->phi() );

    mL1TkMu_trkIsol.push_back( Tkmu->trkIsol() );
    mL1TkMu_trkzVtx.push_back( Tkmu->trkzVtx() );
    mL1TkMu_dR.push_back( Tkmu->dR() );

    mL1TkMu_nTracksMatched.push_back( Tkmu->nTracksMatched() );
    mL1TkMu_trackCurvature.push_back( Tkmu->trackCurvature() );

    mL1TkMu_quality.push_back( Tkmu->quality() );
    mL1TkMu_pattern.push_back( Tkmu->pattern() );
    mL1TkMu_muonDetector.push_back( Tkmu->muonDetector() );

    if (Tkmu->muRef().isNonnull()) {
      auto regionalCandidate = Tkmu->muRef().get();
      mL1TkMu_muRefHwPt.push_back( static_cast<float>(regionalCandidate->hwPt())*0.5 );
      mL1TkMu_muRefHwDXY.push_back( regionalCandidate->hwDXY() ); // 4 bit information, don't know how to decode
      mL1TkMu_muRefHwEta.push_back( static_cast<float>(regionalCandidate->hwEta())*0.010875 );
      mL1TkMu_muRefHwPhi.push_back( static_cast<float>(regionalCandidate->hwPhi())*2.*M_PI/576. );
      mL1TkMu_muRefHwSign.push_back( ( regionalCandidate->hwSign()%2==0 ) ? 1 : -1 );
      mL1TkMu_muRefHwSignValid.push_back( regionalCandidate->hwSignValid() );
      mL1TkMu_muRefHwQual.push_back( regionalCandidate->hwQual() );
    }
    else {
      mL1TkMu_muRefHwPt.push_back( -99999. );
      mL1TkMu_muRefHwDXY.push_back( -99999 );
      mL1TkMu_muRefHwEta.push_back( -99999. );
      mL1TkMu_muRefHwPhi.push_back( -99999. );
      mL1TkMu_muRefHwSign.push_back( -99999 );
      mL1TkMu_muRefHwSignValid.push_back( -99999 );
      mL1TkMu_muRefHwQual.push_back( -99999 );
    }

    auto theTTTrack = Tkmu->trkPtr();
    MuonHLTobjCorrelator::L1TTTrack theTTTrackTmp( *theTTTrack );

    std::map<MuonHLTobjCorrelator::L1TTTrack,unsigned int>::const_iterator where = mTTTrackMap.find(theTTTrackTmp);
    int linkNo = (where==mTTTrackMap.end()) ? -1 : static_cast<int>(mTTTrackMap[theTTTrackTmp]);
    mL1TkMu_TTTpointer.push_back( linkNo );
  }

  edm::Handle<std::vector<l1t::TkPrimaryVertex> > l1TkPrimaryVertex;
  iEvent.getByToken(l1TkPrimaryVertexToken_,l1TkPrimaryVertex);
  for(auto Tkvtx=l1TkPrimaryVertex->begin(); Tkvtx!=l1TkPrimaryVertex->end(); ++Tkvtx)
  {
    // https://github.com/cms-sw/cmssw/blob/09b17fcfb3900782ab78ad6e0c76e1957c94ff71/DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h
    // cout<< "zvertex"<<Tkvtx->zvertex()<<endl;
  }
}
*/

void MuonHLTNtupler::Fill_Muon(const edm::Event &iEvent)
{
  edm::Handle<std::vector<reco::Muon> > h_offlineMuon;
  // edm::Handle< edm::View<reco::Muon> > h_offlineMuon;
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
      // muon_dB_[_nMuon] = mu->dB(); // -- dB is only availabe in pat::Muon -- //
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
  if( pathName.find("HLT_IsoMu")    != std::string::npos ||
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
  if( (filterName.find("sMu") != std::string::npos || filterName.find("SingleMu") != std::string::npos) &&
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
  edm::Handle<reco::GenParticleCollection> h_genParticle;
  iEvent.getByToken(t_genParticle_, h_genParticle);

  int _nGenParticle = 0;
  for( size_t i=0; i< h_genParticle->size(); ++i)
  {
    const reco::GenParticle &parCand = (*h_genParticle)[i];

    if( abs(parCand.pdgId()) == 13 ) // -- only muons -- //
    {
      genParticle_ID_[_nGenParticle]     = parCand.pdgId();
      genParticle_status_[_nGenParticle] = parCand.status();
      // genParticle_mother_[_nGenParticle] = parCand.mother(0)->pdgId();

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
  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  iEvent.getByToken(associatorToken, theAssociator);
  edm::Handle<TrackingParticleCollection> TPCollection;
  iEvent.getByToken(trackingParticleToken, TPCollection);

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

  fill_trackTemplate(iEvent, t_hltIterL3OIMuonTrack_,          theAssociator, TPCollection, tracker,                       hltIterL3OIMuonTrackMap,         TThltIterL3OIMuonTrack);
  fill_trackTemplate(iEvent, t_hltIter0IterL3MuonTrack_,       theAssociator, TPCollection, tracker,       hltIter0IterL3MuonTrackMap,      TThltIter0IterL3MuonTrack);
  fill_trackTemplate(iEvent, t_hltIter2IterL3MuonTrack_,       theAssociator, TPCollection, tracker,                      hltIter2IterL3MuonTrackMap,      TThltIter2IterL3MuonTrack);
  fill_trackTemplate(iEvent, t_hltIter3IterL3MuonTrack_,       theAssociator, TPCollection, tracker,                      hltIter3IterL3MuonTrackMap,      TThltIter3IterL3MuonTrack);
  fill_trackTemplate(iEvent, t_hltIter0IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, hltIter0IterL3FromL1MuonTrackMap,TThltIter0IterL3FromL1MuonTrack);
  fill_trackTemplate(iEvent, t_hltIter2IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker,                hltIter2IterL3FromL1MuonTrackMap,TThltIter2IterL3FromL1MuonTrack);
  fill_trackTemplate(iEvent, t_hltIter3IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker,                hltIter3IterL3FromL1MuonTrackMap,TThltIter3IterL3FromL1MuonTrack);

  // fill_trackTemplate(iEvent, t_hltIterL3OIMuonTrack_,          theAssociator, TPCollection, tracker, mvaHltIterL3OISeedsFromL2Muons_,                       hltIterL3OIMuonTrackMap,         TThltIterL3OIMuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter0IterL3MuonTrack_,       theAssociator, TPCollection, tracker, mvaHltIter0IterL3MuonPixelSeedsFromPixelTracks_,       hltIter0IterL3MuonTrackMap,      TThltIter0IterL3MuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter2IterL3MuonTrack_,       theAssociator, TPCollection, tracker, mvaHltIter2IterL3MuonPixelSeeds_,                      hltIter2IterL3MuonTrackMap,      TThltIter2IterL3MuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter3IterL3MuonTrack_,       theAssociator, TPCollection, tracker, mvaHltIter3IterL3MuonPixelSeeds_,                      hltIter3IterL3MuonTrackMap,      TThltIter3IterL3MuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter0IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, mvaHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_, hltIter0IterL3FromL1MuonTrackMap,TThltIter0IterL3FromL1MuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter2IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, mvaHltIter2IterL3FromL1MuonPixelSeeds_,                hltIter2IterL3FromL1MuonTrackMap,TThltIter2IterL3FromL1MuonTrack);
  // fill_trackTemplate(iEvent, t_hltIter3IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, mvaHltIter3IterL3FromL1MuonPixelSeeds_,                hltIter3IterL3FromL1MuonTrackMap,TThltIter3IterL3FromL1MuonTrack);
}

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
  // pairSeedMvaEstimator pairMvaEstimator,
  std::map<tmpTSOD,unsigned int>& trkMap,
  trkTemplate* TTtrack
) {

  // edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  // bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);

  // edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  // bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

  // edm::Handle<l1t::TkMuonCollection> h_L1TkMu;
  // bool hasL1TkMu = iEvent.getByToken( TkMuonToken_, h_L1TkMu);

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( theToken, trkHandle ) )
  {
    auto recSimColl = theAssociator_->associateRecoToSim(trkHandle,TPCollection_);

    for( unsigned int i = 0; i < trkHandle->size(); i++ )
    {
      TTtrack->fill(trkHandle->at(i));

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

      /*
      if( doMVA && hasL1 && hasL2 && hasL1TkMu ) {
        const TrajectorySeed seed = *(trkHandle->at(i).seedRef());
        GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
        GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

        vector<float> mva = getSeedMva(
          pairMvaEstimator,
          seed,
          global_p,
          global_x,
          h_L1Muon,
          h_L2Muon,
          h_L1TkMu
        );
        TTtrack->fillMva( mva[0], mva[1], mva[2], mva[3] );
      }
      else {
        // cout << "fill_trackTemplate: !(hasL1 && hasL2 && hasL1TkMu)" << endl;
        TTtrack->fillMva( -99999., -99999., -99999., -99999. );
      }
      */
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
  // for( int i=0; i<4; ++i ) {
  //   delete mvaHltIterL3OISeedsFromL2Muons_.at(i).first;
  //   delete mvaHltIterL3OISeedsFromL2Muons_.at(i).second;
  //   delete mvaHltIter0IterL3MuonPixelSeedsFromPixelTracks_.at(i).first;
  //   delete mvaHltIter0IterL3MuonPixelSeedsFromPixelTracks_.at(i).second;
  //   delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).second;
  //   delete mvaHltIter3IterL3MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter3IterL3MuonPixelSeeds_.at(i).second;
  //   delete mvaHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_.at(i).first;
  //   delete mvaHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_.at(i).second;
  //   delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).second;
  //   delete mvaHltIter3IterL3FromL1MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter3IterL3FromL1MuonPixelSeeds_.at(i).second;
  // }
}

void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
