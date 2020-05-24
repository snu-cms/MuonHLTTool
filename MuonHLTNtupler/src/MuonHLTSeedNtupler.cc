// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTSeedNtupler.h"

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
#include "DataFormats/Math/interface/deltaPhi.h"
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


MuonHLTSeedNtupler::MuonHLTSeedNtupler(const edm::ParameterSet& iConfig):
// trackerHitAssociatorConfig_(iConfig, consumesCollector()),
associatorToken(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("associator"))),
trackingParticleToken(consumes<TrackingParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingParticle"))),

t_offlineVertex_     ( consumes< reco::VertexCollection >                 (iConfig.getUntrackedParameter<edm::InputTag>("offlineVertex"     )) ),
t_PUSummaryInfo_     ( consumes< std::vector<PileupSummaryInfo> >         (iConfig.getUntrackedParameter<edm::InputTag>("PUSummaryInfo"     )) ),

t_L1Muon_            ( consumes< l1t::MuonBxCollection  >                 (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon"            )) ),
t_L2Muon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2Muon"            )) ),

t_L1TkMuon_          ( consumes< l1t::TkMuonCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("L1TkMuon"))),
t_L1TkPrimaryVertex_ ( consumes< l1t::TkPrimaryVertexCollection >         (iConfig.getUntrackedParameter<edm::InputTag>("L1TkPrimaryVertex"))),

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

t_genParticle_       ( consumes< reco::GenParticleCollection >            (iConfig.getUntrackedParameter<edm::InputTag>("genParticle"       )) )

{}

void MuonHLTSeedNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  testRun(iEvent, iSetup, t_hltIter2IterL3FromL1MuonPixelSeeds_);

  // Init();

  // -- fill each object
  // Fill_Event(iEvent);
  // Fill_IterL3TT(iEvent);
  // Fill_Seed(iEvent, iSetup);
  // NTEvent_->Fill();
}

void MuonHLTSeedNtupler::beginJob()
{
  edm::Service<TFileService> fs;

  NTEvent_    = fs->make<TTree>("NTEvent","NTEvent");  

  NThltIterL3OI_    = fs->make<TTree>("NThltIterL3OI","NThltIterL3OI");

  NThltIter0_       = fs->make<TTree>("NThltIter0","NThltIter0");
  NThltIter2_       = fs->make<TTree>("NThltIter2","NThltIter2");
  NThltIter3_       = fs->make<TTree>("NThltIter3","NThltIter3");

  NThltIter0FromL1_ = fs->make<TTree>("NThltIter0FromL1","NThltIter0FromL1");
  NThltIter2FromL1_ = fs->make<TTree>("NThltIter2FromL1","NThltIter2FromL1");
  NThltIter3FromL1_ = fs->make<TTree>("NThltIter3FromL1","NThltIter3FromL1");

  Make_Branch();

  // HERE
  t = fs->make<TTree>("t","t");
  t->Branch("l1tt_tsos_x",  &l1tt_tsos_x_, "l1tt_tsos_x/F");
  t->Branch("l1tt_tsos_y",  &l1tt_tsos_y_, "l1tt_tsos_y/F");
  t->Branch("l1tt_tsos_z",  &l1tt_tsos_z_, "l1tt_tsos_z/F");
  t->Branch("hit_x",  &hit_x_, "hit_x/F");
  t->Branch("hit_y",  &hit_y_, "hit_y/F");
  t->Branch("hit_z",  &hit_z_, "hit_z/F");
}

void MuonHLTSeedNtupler::Init()
{
  runNum_       = -999;
  lumiBlockNum_ = -999;
  eventNum_     = 0;
  nVertex_      = -999;
  truePU_       = -999;
  nhltIterL3OI_ = 0;
  nhltIter0_ = 0;
  nhltIter2_ = 0;
  nhltIter3_ = 0;
  nhltIter0FromL1_ = 0;
  nhltIter2FromL1_ = 0;
  nhltIter3FromL1_ = 0;
  hltIterL3OIMuonTrackMap.clear();
  hltIter0IterL3MuonTrackMap.clear();
  hltIter2IterL3MuonTrackMap.clear();
  hltIter3IterL3MuonTrackMap.clear();
  hltIter0IterL3FromL1MuonTrackMap.clear();
  hltIter2IterL3FromL1MuonTrackMap.clear();
  hltIter3IterL3FromL1MuonTrackMap.clear();

  ST->clear();
  // SThltIterL3OISeedsFromL2Muons->clear();
  // SThltIter0IterL3MuonPixelSeedsFromPixelTracks->clear();
  // SThltIter2IterL3MuonPixelSeeds->clear();
  // SThltIter3IterL3MuonPixelSeeds->clear();
  // SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->clear();
  // SThltIter2IterL3FromL1MuonPixelSeeds->clear();
  // SThltIter3IterL3FromL1MuonPixelSeeds->clear();

  TThltIterL3OIMuonTrack->clear();
  TThltIter0IterL3MuonTrack->clear();
  TThltIter2IterL3MuonTrack->clear();
  TThltIter3IterL3MuonTrack->clear();
  TThltIter0IterL3FromL1MuonTrack->clear();
  TThltIter2IterL3FromL1MuonTrack->clear();
  TThltIter3IterL3FromL1MuonTrack->clear();
}

void MuonHLTSeedNtupler::Make_Branch()
{
  NTEvent_->Branch("runNum",&runNum_,"runNum/I");
  NTEvent_->Branch("lumiBlockNum",&lumiBlockNum_,"lumiBlockNum/I");
  NTEvent_->Branch("eventNum",&eventNum_,"eventNum/l"); // -- unsigned long long -- //
  NTEvent_->Branch("nVertex", &nVertex_, "nVertex/I");
  NTEvent_->Branch("truePU", &truePU_, "truePU/I");
  NTEvent_->Branch("nhltIterL3OI",  &nhltIterL3OI_, "nhltIterL3OI/I");
  NTEvent_->Branch("nhltIter0",  &nhltIter0_, "nhltIter0/I");
  NTEvent_->Branch("nhltIter2",  &nhltIter2_, "nhltIter2/I");
  NTEvent_->Branch("nhltIter3",  &nhltIter3_, "nhltIter3/I");
  NTEvent_->Branch("nhltIter0FromL1",  &nhltIter0FromL1_, "nhltIter0FromL1/I");
  NTEvent_->Branch("nhltIter2FromL1",  &nhltIter2FromL1_, "nhltIter2FromL1/I");
  NTEvent_->Branch("nhltIter3FromL1",  &nhltIter3FromL1_, "nhltIter3FromL1/I");
  ST->setBranch(NThltIterL3OI_);
  ST->setBranch(NThltIter0_);
  ST->setBranch(NThltIter2_);
  ST->setBranch(NThltIter3_);
  ST->setBranch(NThltIter0FromL1_);
  ST->setBranch(NThltIter2FromL1_);
  ST->setBranch(NThltIter3FromL1_);

  // SThltIterL3OISeedsFromL2Muons->setBranch(NThltIterL3OI_);
  // SThltIter0IterL3MuonPixelSeedsFromPixelTracks->setBranch(NThltIter0_);
  // SThltIter2IterL3MuonPixelSeeds->setBranch(NThltIter2_);
  // SThltIter3IterL3MuonPixelSeeds->setBranch(NThltIter3_);
  // SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->setBranch(NThltIter0FromL1_);
  // SThltIter2IterL3FromL1MuonPixelSeeds->setBranch(NThltIter2FromL1_);
  // SThltIter3IterL3FromL1MuonPixelSeeds->setBranch(NThltIter3FromL1_);
}

void MuonHLTSeedNtupler::Fill_Event(const edm::Event &iEvent)
{
  // -- basic info.
  runNum_       = iEvent.id().run();
  lumiBlockNum_ = iEvent.id().luminosityBlock();
  eventNum_     = iEvent.id().event();

  // -- pile up
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

  // -- vertex
  edm::Handle<reco::VertexCollection> h_offlineVertex;
  if( iEvent.getByToken(t_offlineVertex_, h_offlineVertex) )
  {
    int nGoodVtx = 0;
    for(reco::VertexCollection::const_iterator it = h_offlineVertex->begin(); it != h_offlineVertex->end(); ++it)
      if( it->isValid() ) nGoodVtx++;

    nVertex_ = nGoodVtx;
  }
}

void MuonHLTSeedNtupler::Fill_IterL3TT(const edm::Event &iEvent)
{
  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  edm::Handle<TrackingParticleCollection> TPCollection;

  if( iEvent.getByToken(associatorToken, theAssociator) && iEvent.getByToken(trackingParticleToken, TPCollection) ) {
    fill_trackTemplate(iEvent,t_hltIterL3OIMuonTrack_,theAssociator,TPCollection,hltIterL3OIMuonTrackMap,TThltIterL3OIMuonTrack);
    fill_trackTemplate(iEvent,t_hltIter0IterL3MuonTrack_,theAssociator,TPCollection,hltIter0IterL3MuonTrackMap,TThltIter0IterL3MuonTrack);
    fill_trackTemplate(iEvent,t_hltIter2IterL3MuonTrack_,theAssociator,TPCollection,hltIter2IterL3MuonTrackMap,TThltIter2IterL3MuonTrack);
    fill_trackTemplate(iEvent,t_hltIter3IterL3MuonTrack_,theAssociator,TPCollection,hltIter3IterL3MuonTrackMap,TThltIter3IterL3MuonTrack);
    fill_trackTemplate(iEvent,t_hltIter0IterL3FromL1MuonTrack_,theAssociator,TPCollection,hltIter0IterL3FromL1MuonTrackMap,TThltIter0IterL3FromL1MuonTrack);
    fill_trackTemplate(iEvent,t_hltIter2IterL3FromL1MuonTrack_,theAssociator,TPCollection,hltIter2IterL3FromL1MuonTrackMap,TThltIter2IterL3FromL1MuonTrack);
    fill_trackTemplate(iEvent,t_hltIter3IterL3FromL1MuonTrack_,theAssociator,TPCollection,hltIter3IterL3FromL1MuonTrackMap,TThltIter3IterL3FromL1MuonTrack);
  }
}

void MuonHLTSeedNtupler::Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // TrackerHitAssociator associate(iEvent, trackerHitAssociatorConfig_);
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

  fill_seedTemplate(iEvent, t_hltIterL3OISeedsFromL2Muons_, tracker, hltIterL3OIMuonTrackMap, TThltIterL3OIMuonTrack, NThltIterL3OI_, nhltIterL3OI_ );
  fill_seedTemplate(iEvent, t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_, tracker, hltIter0IterL3MuonTrackMap, TThltIter0IterL3MuonTrack, NThltIter0_, nhltIter0_ );
  fill_seedTemplate(iEvent, t_hltIter2IterL3MuonPixelSeeds_, tracker, hltIter2IterL3MuonTrackMap, TThltIter2IterL3MuonTrack, NThltIter2_, nhltIter2_ );
  fill_seedTemplate(iEvent, t_hltIter3IterL3MuonPixelSeeds_, tracker, hltIter3IterL3MuonTrackMap, TThltIter3IterL3MuonTrack, NThltIter3_, nhltIter3_ );
  fill_seedTemplate(iEvent, t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_, tracker, hltIter0IterL3FromL1MuonTrackMap, TThltIter0IterL3FromL1MuonTrack, NThltIter0FromL1_, nhltIter0FromL1_ );
  fill_seedTemplate(iEvent, t_hltIter2IterL3FromL1MuonPixelSeeds_, tracker, hltIter2IterL3FromL1MuonTrackMap, TThltIter2IterL3FromL1MuonTrack, NThltIter2FromL1_, nhltIter2FromL1_ );
  fill_seedTemplate(iEvent, t_hltIter3IterL3FromL1MuonPixelSeeds_, tracker, hltIter3IterL3FromL1MuonTrackMap, TThltIter3IterL3FromL1MuonTrack, NThltIter3FromL1_, nhltIter3FromL1_ );
}

void MuonHLTSeedNtupler::fill_trackTemplate(const edm::Event &iEvent, edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
  edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_, edm::Handle<TrackingParticleCollection>& TPCollection_,
  std::map<tmpTSOD,unsigned int>& trkMap, trkTemplate* TTtrack) {

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( theToken, trkHandle ) )
  {
    auto recSimColl = theAssociator_->associateRecoToSim(trkHandle,TPCollection_);

    for( unsigned int i = 0; i < trkHandle->size(); i++ )
    {
      TTtrack->fill(trkHandle->at(i));

      int linkNo = -1;
      // for (unsigned int idxL3passed = 0; idxL3passed < iterL3IDpassed.size(); idxL3passed++) {
      //   if ( iterL3IDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNo = idxL3passed;
      // }
      TTtrack->linkIterL3(linkNo);

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
      } else {
        TTtrack->fillDummyTP();
        TTtrack->fillBestTPsharedFrac(-99999.);
        TTtrack->fillmatchedTPsize(0);
      }
    }
  }
}

void MuonHLTSeedNtupler::fill_seedTemplate(
  const edm::Event &iEvent, edm::EDGetTokenT<TrajectorySeedCollection>& theToken,
  edm::ESHandle<TrackerGeometry>& tracker, std::map<tmpTSOD,unsigned int>& trkMap, trkTemplate* TTtrack, TTree* NT, int &nSeed ) {

  edm::Handle<reco::GenParticleCollection> h_genParticle;
  bool hasGen = iEvent.getByToken(t_genParticle_, h_genParticle);

  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  bool hasL1 = iEvent.getByToken(t_L1Muon_, h_L1Muon);

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

  edm::Handle<l1t::TkMuonCollection> h_L1TkMu;
  bool hasL1TkMu = iEvent.getByToken(t_L1TkMuon_, h_L1TkMu);

  // edm::Handle<l1t::TkPrimaryVertexCollection> h_L1TkPrimaryVertex;
  // bool hasL1TkVtx = iEvent.getByToken(t_L1TkPrimaryVertex_, h_L1TkPrimaryVertex);

  edm::Handle< TrajectorySeedCollection > seedHandle;
  if( iEvent.getByToken( theToken, seedHandle) )
  {
    nSeed = seedHandle->size();
    for( auto i=0U; i<seedHandle->size(); ++i )
    {
      const auto& seed(seedHandle->at(i));

      tmpTSOD seedTsod(seed.startingState());
      ST->clear();

      ST->fill_PU(
        truePU_
      );

      // -- Track association
      ST->fill(seed, tracker);
      std::map<tmpTSOD,unsigned int>::const_iterator where = trkMap.find(seedTsod);
      int idxtmpL3 = (where==trkMap.end()) ? -1 : trkMap[seedTsod];
      ST->fill_TP(TTtrack, idxtmpL3 );

      GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
      GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

      // -- GenParticle (muon) tag -- //
      if( hasGen )
      {
        float gen_pt = -99999.;
        float gen_eta = -99999.;
        float gen_phi = -99999.;
        float dR_GenSeed = 99999.;
        for(auto genp = h_genParticle->begin(); genp != h_genParticle->end(); genp++)
        {
          if( fabs(genp->pdgId()) ==  13 && genp->status()==1 )
          {
            GlobalVector vec_seed_vtx( global_x.x() - genp->vx(), global_x.y() - genp->vy(), global_x.z() - genp->vz() );
            if( reco::deltaR( *genp, vec_seed_vtx ) < dR_GenSeed ) {
              dR_GenSeed = reco::deltaR( *genp, vec_seed_vtx );
              gen_pt = genp->pt();
              gen_eta = genp->eta();
              gen_phi = genp->phi();
            }
          }
        }
        ST->fill_Genvars(
          gen_pt,
          gen_eta,
          gen_phi
        );
      }

      // -- L1, L2 association
      if( hasL1 ) {
        float dR_minDRL1SeedP = 99999.;
        float dPhi_minDRL1SeedP = 99999.;
        float dR_minDPhiL1SeedX = 99999.;
        float dPhi_minDPhiL1SeedX = 99999.;
        float dR_minDRL1SeedP_AtVtx = 99999.;
        float dPhi_minDRL1SeedP_AtVtx = 99999.;
        float dR_minDPhiL1SeedX_AtVtx = 99999.;
        float dPhi_minDPhiL1SeedX_AtVtx = 99999.;
        for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
        {
          if(ibx != 0) continue; // -- only take when ibx == 0 -- //
          for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
          {
            l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

            if(ref_L1Mu->hwQual() < 7)
              continue;

            float dR_L1SeedP   = reco::deltaR( *ref_L1Mu, global_p);
            float dPhi_L1SeedP = reco::deltaPhi( ref_L1Mu->phi(), global_p.phi());
            float dR_L1SeedX   = reco::deltaR( *ref_L1Mu, global_x);
            float dPhi_L1SeedX = reco::deltaPhi( ref_L1Mu->phi(), global_x.phi());

            if( dR_L1SeedP < dR_minDRL1SeedP ) {
              dR_minDRL1SeedP = dR_L1SeedP;
              dPhi_minDRL1SeedP = dPhi_L1SeedP;
            }
            if( fabs(dPhi_L1SeedX) < fabs(dPhi_minDPhiL1SeedX) ) {
              dR_minDPhiL1SeedX = dR_L1SeedX;
              dPhi_minDPhiL1SeedX = dPhi_L1SeedX;
            }

            float dR_L1SeedP_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_p.eta(), global_p.phi());
            float dPhi_L1SeedP_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_p.phi());
            float dR_L1SeedX_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_x.eta(), global_x.phi());
            float dPhi_L1SeedX_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_x.phi());

            if( dR_L1SeedP_AtVtx < dR_minDRL1SeedP_AtVtx ) {
              dR_minDRL1SeedP_AtVtx = dR_L1SeedP_AtVtx;
              dPhi_minDRL1SeedP_AtVtx = dPhi_L1SeedP_AtVtx;
            }
            if( fabs(dPhi_L1SeedX_AtVtx) < fabs(dPhi_minDPhiL1SeedX_AtVtx) ) {
              dR_minDPhiL1SeedX_AtVtx = dR_L1SeedX_AtVtx;
              dPhi_minDPhiL1SeedX_AtVtx = dPhi_L1SeedX_AtVtx;
            }
          }
        }
        ST->fill_L1vars(
          dR_minDRL1SeedP,         dPhi_minDRL1SeedP,
          dR_minDPhiL1SeedX,       dPhi_minDPhiL1SeedX,
          dR_minDRL1SeedP_AtVtx,   dPhi_minDRL1SeedP_AtVtx,
          dR_minDPhiL1SeedX_AtVtx, dPhi_minDPhiL1SeedX_AtVtx
        );
      }

      if( hasL2 && h_L2Muon->size() > 0 ) {
        float dR_minDRL2SeedP = 99999.;
        float dPhi_minDRL2SeedP = 99999.;
        float dR_minDPhiL2SeedX = 99999.;
        float dPhi_minDPhiL2SeedX = 99999.;
        for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
        {
          reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, i_L2);

          float dR_L2SeedP   = reco::deltaR( *ref_L2Mu, global_p);
          float dPhi_L2SeedP = reco::deltaPhi( ref_L2Mu->phi(), global_p.phi());
          float dR_L2SeedX   = reco::deltaR( *ref_L2Mu, global_x);
          float dPhi_L2SeedX = reco::deltaPhi( ref_L2Mu->phi(), global_x.phi());

          if( dR_L2SeedP < dR_minDRL2SeedP ) {
            dR_minDRL2SeedP = dR_L2SeedP;
            dPhi_minDRL2SeedP = dPhi_L2SeedP;
          }
          if( fabs(dPhi_L2SeedX) < fabs(dPhi_minDPhiL2SeedX) ) {
            dR_minDPhiL2SeedX = dR_L2SeedX;
            dPhi_minDPhiL2SeedX = dPhi_L2SeedX;
          }
        }

        ST->fill_L2vars(
          dR_minDRL2SeedP,         dPhi_minDRL2SeedP,
          dR_minDPhiL2SeedX,       dPhi_minDPhiL2SeedX
        );
      }

      // -- L1TkMu association
      if( hasL1TkMu && h_L1TkMu->size() > 0 ) {
        float dR_L1TkMuSeedP = 99999.;
        float dPhi_L1TkMuSeedP = 99999.;
        for(auto L1TkMu=h_L1TkMu->begin(); L1TkMu!=h_L1TkMu->end(); ++L1TkMu)
        {
          auto TkRef = L1TkMu->trkPtr();
          float dR_L1TkMuSeedP_tmp   = reco::deltaR( *TkRef, global_p);
          float dPhi_L1TkMuSeedP_tmp = reco::deltaPhi( TkRef->phi(), global_p.phi());
          if( dR_L1TkMuSeedP_tmp < dR_L1TkMuSeedP ) {
            dR_L1TkMuSeedP   = dR_L1TkMuSeedP_tmp;
            dPhi_L1TkMuSeedP = dPhi_L1TkMuSeedP_tmp;
          }
        }
        ST->fill_L1TkMuvars(
          dR_L1TkMuSeedP,
          dPhi_L1TkMuSeedP
        );
      }

      // --L1TkMu, L1TkVertex association
      // if( hasL1TkVtx ) {
      //   for(auto L1TkVtx=h_L1TkPrimaryVertex->begin(); L1TkVtx!=h_L1TkPrimaryVertex->end(); ++L1TkVtx)
      //   {
      //   }
      // }

      ST->fill_ntuple(NT);
    } // -- end of seed iteration
  } // -- if getByToken is valid
}

void MuonHLTSeedNtupler::testRun(
  const edm::Event &iEvent, const edm::EventSetup& iSetup,
  edm::EDGetTokenT<TrajectorySeedCollection>& theToken
) {
  edm::ESHandle<MagneticField> magfieldH;
  iSetup.get<IdealMagneticFieldRecord>().get(magfieldH);

  edm::ESHandle<Propagator> propagatorAlongH;
  iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterialParabolicMf", propagatorAlongH);
  std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);

  edm::ESHandle<TrackerGeometry> trkgeom;
  iSetup.get<TrackerDigiGeometryRecord>().get(trkgeom);

  edm::ESHandle<GeometricDet> geomDet;
  iSetup.get<IdealGeometryRecord>().get(geomDet);

  edm::ESHandle<TrackerTopology> trkTopo;
  iSetup.get<TrackerTopologyRcd>().get(trkTopo);

  GeometricSearchTrackerBuilder builder;
  GeometricSearchTracker* geomTracker = builder.build(&(*geomDet), &(*trkgeom), &(*trkTopo));

  edm::Handle<l1t::TkMuonCollection> h_L1TkMu;
  bool hasL1TkMu = iEvent.getByToken(t_L1TkMuon_, h_L1TkMu);

  edm::Handle< TrajectorySeedCollection > seedHandle;
  if( hasL1TkMu && iEvent.getByToken( theToken, seedHandle) ) {
    for( auto i=0U; i<seedHandle->size(); ++i ) {
      // HERE
      l1tt_tsos_x_ = -99999.;
      l1tt_tsos_y_ = -99999.;
      l1tt_tsos_z_ = -99999.;
      hit_x_ = -99999.;
      hit_y_ = -99999.;
      hit_z_ = -99999.;

      const auto& seed(seedHandle->at(i));

      vector< pair<LayerHit, LayerTSOS> > hitTsosPairs = getHitTsosPairs(
        seed,
        h_L1TkMu,
        magfieldH,
        *(propagatorAlong.get()),
        geomTracker
      );

      // HERE
      bool found = (hitTsosPairs.size()>0);
      cout << i << ": Matched found? " << found << endl;
      for( auto it=hitTsosPairs.begin(); it!=hitTsosPairs.end(); ++it ) {
        // auto lhit  = it->first.first;
        auto hit   = it->first.second;
        // auto ltsos = it->second.first;
        auto tsos  = it->second.second;

        float deltaz = hit->globalPosition().z() - tsos.globalPosition().z();
        float deltaxy = ( hit->globalPosition() - tsos.globalPosition() ).perp();
        float dr = ( hit->globalPosition() - tsos.globalPosition() ).mag();
        float dR = reco::deltaR( tsos.globalPosition(), hit->globalPosition() );

        cout << "\t the hit:  " << hit->globalPosition() << "\t" << hit->globalPosition().perp() << endl;
        cout << "\t the tsos: " << tsos.globalPosition() << "\t" << tsos.globalPosition().perp() << endl;
        cout << "\t dr= " << dr << "\t delta z: " << deltaz << "\t delta r: " << deltaxy << "\t dR= " << dR << endl;

        hit_x_ = hit->globalPosition().x();
        hit_y_ = hit->globalPosition().y();
        hit_z_ = hit->globalPosition().z();
        l1tt_tsos_x_ = tsos.globalPosition().x();
        l1tt_tsos_y_ = tsos.globalPosition().y();
        l1tt_tsos_z_ = tsos.globalPosition().z();
        t->Fill();
      }

    }
  }

  return;
}

vector< LayerTSOS > MuonHLTSeedNtupler::getTsosOnPixels(
  l1t::TkMuon L1TkMu,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker
) {
  vector< LayerTSOS > v_tsos = {};

  std::vector<BarrelDetLayer const*>  const&  bpix = geomTracker->pixelBarrelLayers();
  std::vector<ForwardDetLayer const*> const& nfpix = geomTracker->negPixelForwardLayers();
  std::vector<ForwardDetLayer const*> const& pfpix = geomTracker->posPixelForwardLayers();

  // -- L1TkMu selection
  // if( L1TkMu->muRef().isNull() )  continue;
  if( L1TkMu.trkPtr().isNull() )  return v_tsos;
  // FIXME this is random choice

  auto l1tk = *(L1TkMu.trkPtr());
  int chargeTk = l1tk.rInv() > 0. ? 1 : -1;
  GlobalPoint  gpos = l1tk.POCA();
  GlobalVector gmom = l1tk.momentum();

  FreeTrajectoryState fts = FreeTrajectoryState( gpos, gmom, chargeTk, magfieldH.product() );

  for(auto it = bpix.begin(); it != bpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto z0 = std::abs(tsos.globalPosition().z() - (**it).specificSurface().position().z());
    auto deltaZ = 0.5f * (**it).surface().bounds().length();
    deltaZ += 0.5f * (**it).surface().bounds().thickness() * std::abs(tsos.globalDirection().z()) / tsos.globalDirection().perp();
    bool is_compatible  = (z0 < deltaZ);

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }
  for(auto it = nfpix.begin(); it != nfpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto r2 = tsos.localPosition().perp2();
    float deltaR = 0.5f * (**it).surface().bounds().thickness() * tsos.localDirection().perp() / std::abs(tsos.localDirection().z());
    auto ri2 = std::max((**it).specificSurface().innerRadius() - deltaR, 0.f);
    ri2 *= ri2;
    auto ro2 = (**it).specificSurface().outerRadius() + deltaR;
    ro2 *= ro2;
    bool is_compatible  = ((r2 > ri2) && (r2 < ro2));

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }
  for(auto it = pfpix.begin(); it != pfpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto r2 = tsos.localPosition().perp2();
    float deltaR = 0.5f * (**it).surface().bounds().thickness() * tsos.localDirection().perp() / std::abs(tsos.localDirection().z());
    auto ri2 = std::max((**it).specificSurface().innerRadius() - deltaR, 0.f);
    ri2 *= ri2;
    auto ro2 = (**it).specificSurface().outerRadius() + deltaR;
    ro2 *= ro2;
    bool is_compatible  = ((r2 > ri2) && (r2 < ro2));

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }

  return v_tsos;
}

// -- hit, TSOS pairs for each L1TkMu
vector< pair<LayerHit, LayerTSOS> > MuonHLTSeedNtupler::getHitTsosPairs(
  TrajectorySeed seed,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker
) {
  vector< pair<LayerHit, LayerTSOS> > out = {};

  // FIXME this is random choice
  float av_dr_min = 20.;

  // -- loop on L1TkMu
  for(auto L1TkMu=h_L1TkMu->begin(); L1TkMu!=h_L1TkMu->end(); ++L1TkMu) {

    vector< LayerTSOS > v_tsos = getTsosOnPixels(
      *L1TkMu,
      magfieldH,
      propagatorAlong,
      geomTracker
    );

    // -- loop on recHits
    vector<int> v_tsos_skip( v_tsos.size(), 0 );
    vector< pair<LayerHit, LayerTSOS> > hitTsosPair = {};
    int ihit = 0;
    for( auto hit = seed.recHits().first; hit!=seed.recHits().second; ++hit ) {
      // -- look for closest tsos by absolute distance
      // FIXME this is random choice
      int the_tsos = -99999;
      float dr_min = 20.;
      for( auto i=0U; i<v_tsos.size(); ++i ) {
        if( v_tsos_skip.at(i) )  continue;
        float dr = ( v_tsos.at(i).second.globalPosition() - hit->globalPosition() ).mag();
        if( dr < dr_min ) {
          dr_min = dr;
          the_tsos = i;
          v_tsos_skip.at(i) = 1;
        }
      }

      if( the_tsos > -1 ) {
        const DetLayer* thelayer =  geomTracker->idToLayer( hit->geographicalId() );
        hitTsosPair.push_back( make_pair( make_pair( thelayer, &*hit), v_tsos.at(the_tsos) ) );
      }

      ihit++;
    } // loop on recHits

    // -- find tsos for all recHits?
    // FIXME this is random choice
    if( (int)hitTsosPair.size() == ihit ) {
      float av_dr = 0.;
      for( auto it=hitTsosPair.begin(); it!=hitTsosPair.end(); ++it ) {
        auto hit  = it->first.second;
        auto tsos = it->second.second;
        av_dr += ( hit->globalPosition() - tsos.globalPosition() ).mag();
      }
      av_dr = av_dr > 0 ? av_dr / (float)hitTsosPair.size() : 1.e6;

      if( av_dr < av_dr_min ) {
        av_dr_min = av_dr;
        out = hitTsosPair;
      }
    }

  }  // loop on L1TkMu

  return out;
}

void MuonHLTSeedNtupler::endJob() {}
void MuonHLTSeedNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void MuonHLTSeedNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTSeedNtupler);
