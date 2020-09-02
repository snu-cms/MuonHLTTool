#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

class GenMuAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenMuAnalyzer(const edm::ParameterSet&);
  virtual ~GenMuAnalyzer() {};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT< reco::GenParticleCollection >     t_genParticle;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > t_L1TT;
  edm::EDGetTokenT<l1t::TkMuonCollection> t_L1TkMuon;
  double _pt_min;

  TH1F *h_gen_pt;
  TH1F *h_gen_eta;
  TH1F *h_gen_phi_abseta1p2;
  TH1F *h_gen_matL1TkMuon_pt;
  TH1F *h_gen_matL1TkMuon_eta;
  TH1F *h_gen_matL1TkMuon_phi_abseta1p2;
  TH1F *h_gen_matL1TT_pt;
  TH1F *h_gen_matL1TT_eta;
  TH1F *h_gen_matL1TT_phi_abseta1p2;
};


GenMuAnalyzer::GenMuAnalyzer(const edm::ParameterSet& iConfig)
  : t_genParticle( consumes< reco::GenParticleCollection >( iConfig.getParameter<edm::InputTag>("genParticle_src") ) ),
    t_L1TT(        consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( iConfig.getParameter<edm::InputTag>("L1TT_src") ) ),
    t_L1TkMuon(    consumes< l1t::TkMuonCollection >( iConfig.getParameter<edm::InputTag>("L1TkMuon_src") ) ),
    _pt_min(       iConfig.getParameter< double >("pt_min"))
{
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  h_gen_pt                        = fs->make<TH1F>("h_gen_pt",  "", 1000, 0, 1000);
  h_gen_eta                       = fs->make<TH1F>("h_gen_eta", "", 60, -3, 3);
  h_gen_phi_abseta1p2             = fs->make<TH1F>("h_gen_phi_abseta1p2", "", 60, -3.2, 3.2);

  h_gen_matL1TkMuon_pt            = fs->make<TH1F>("h_gen_matL1TkMuon_pt",  "", 1000, 0, 1000);
  h_gen_matL1TkMuon_eta           = fs->make<TH1F>("h_gen_matL1TkMuon_eta", "", 60, -3, 3);
  h_gen_matL1TkMuon_phi_abseta1p2 = fs->make<TH1F>("h_gen_matL1TkMuon_phi_abseta1p2", "", 60, -3.2, 3.2);

  h_gen_matL1TT_pt                = fs->make<TH1F>("h_gen_matL1TT_pt",  "", 1000, 0, 1000);
  h_gen_matL1TT_eta               = fs->make<TH1F>("h_gen_matL1TT_eta", "", 60, -3, 3);
  h_gen_matL1TT_phi_abseta1p2     = fs->make<TH1F>("h_gen_matL1TT_phi_abseta1p2", "", 60, -3.2, 3.2);
}

void GenMuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  double dRcone = 0.3;

  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > h_L1TT;
  bool hasL1TT = iEvent.getByToken(t_L1TT, h_L1TT);

  edm::Handle<l1t::TkMuonCollection> h_L1TkMuon;
  bool hasL1TkMuon = iEvent.getByToken(t_L1TkMuon, h_L1TkMuon);

  edm::Handle<reco::GenParticleCollection> h_genParticle;
  if( hasL1TT && hasL1TkMuon && iEvent.getByToken(t_genParticle,h_genParticle) ) {
    reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
    for (; genp != h_genParticle->end(); genp++) {

      if( fabs(genp->pdgId()) != 13 )  continue;
      if( !genp->isPromptFinalState() )  continue;
      if( !genp->fromHardProcessFinalState() )  continue;
      if( fabs(genp->eta()) > 2.4 )  continue;

      h_gen_pt->Fill( genp->pt() );
      if( genp->pt() > _pt_min )
        h_gen_eta->Fill( genp->eta() );
      if( genp->pt() > _pt_min && fabs(genp->eta()) > 1.2 )
        h_gen_phi_abseta1p2->Fill( genp->phi() );

      bool L1TTmatched = false;
      for(auto itL1TT = h_L1TT->begin(); itL1TT != h_L1TT->end(); itL1TT++) {
        if( reco::deltaR( itL1TT->momentum(), *genp ) < dRcone ) {
          L1TTmatched = true;
          break;
        }
      }
      if( L1TTmatched ) {
        h_gen_matL1TT_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )
          h_gen_matL1TT_eta->Fill( genp->eta() );
        if( genp->pt() > _pt_min && fabs(genp->eta()) > 1.2 )
          h_gen_matL1TT_phi_abseta1p2->Fill( genp->phi() );
      }

      bool L1TkMuonmatched = false;
      for(auto itL1TkMuon = h_L1TkMuon->begin(); itL1TkMuon != h_L1TkMuon->end(); itL1TkMuon++) {
        if( reco::deltaR( itL1TkMuon->momentum(), *genp ) < dRcone ) {
          L1TkMuonmatched = true;
          break;
        }
      }
      if( L1TkMuonmatched ) {
        h_gen_matL1TkMuon_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )
          h_gen_matL1TkMuon_eta->Fill( genp->eta() );
        if( genp->pt() > _pt_min && fabs(genp->eta()) > 1.2 )
          h_gen_matL1TkMuon_phi_abseta1p2->Fill( genp->phi() );
      }

    }
  }

}

DEFINE_FWK_MODULE(GenMuAnalyzer);
