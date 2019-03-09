// -- To investigate the increase of rates in MET triggers
// -- e.g. HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v8

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/Common/interface/TriggerNames.h"

// -- data formats
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <vector>
#include <iostream>

using namespace std;
using namespace edm;

class MHPrintBTagMuTriggerInfo : public edm::EDAnalyzer {
public:
  MHPrintBTagMuTriggerInfo(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

  void printRecoTrackInfo(const edm::Handle < vector<reco::Track> >&, const std::string&);
  void printMuonInfo     (const edm::Handle<  vector<reco::Muon>  >&, const std::string&);
  bool isGoodTrack       (const reco::Track&);

private:
	edm::EDGetTokenT< edm::TriggerResults          > token_trigResults;  
	edm::EDGetTokenT< trigger::TriggerEvent        > token_trigEvent;

	// -- muon track
	edm::EDGetTokenT< vector<reco::Track>          > token_muonMergedTrack;
	edm::EDGetTokenT< vector<reco::Track>          > token_hltIterL3MuonTracks;

	edm::EDGetTokenT< vector<reco::CaloJet>        > caloJet_;

	edm::EDGetTokenT< edm::View<reco::Track >      > softMuTrackPt5_;

	// -- beam spot
	edm::EDGetTokenT< reco::BeamSpot               > beamSpot_;
	reco::Track::Point vertex_;

	// -- muons
	// edm::EDGetTokenT< vector<reco::Muon>           > token_L3MuonsNoID;
	// edm::EDGetTokenT< vector<reco::Muon>           > token_L3Muons;
};


MHPrintBTagMuTriggerInfo::MHPrintBTagMuTriggerInfo(const edm::ParameterSet& iConfig):
token_trigResults           ( consumes< edm::TriggerResults          > ( iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")) ),
token_trigEvent             ( consumes< trigger::TriggerEvent        > ( iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"))   ),
token_muonMergedTrack       ( consumes< vector<reco::Track>          > ( edm::InputTag("hltIterL3MuonAndMuonFromL1Merged::MYHLT"))       ),
token_hltIterL3MuonTracks   ( consumes< vector<reco::Track>          > ( edm::InputTag("hltIterL3MuonTracks::MYHLT"))                    ),
caloJet_                    ( consumes< vector<reco::CaloJet>        > ( edm::InputTag("hltBSoftMuonDiJet40L1FastJetL25Jets::MYHLT"))    ),
softMuTrackPt5_             ( consumes< edm::View<reco::Track >      > ( edm::InputTag("hltBSoftMuonMu5L3::MYHLT"))                      ),
beamSpot_                   ( consumes< reco::BeamSpot               > ( edm::InputTag("hltOnlineBeamSpot::MYHLT"))                      )
{

}

void MHPrintBTagMuTriggerInfo::analyze( const edm::Event& iEvent, const edm::EventSetup& iEventSetup) {

	// -- print run:lumi:event number
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event() << endl;
	cout << endl;

	// -- trigger info.
	edm::Handle<edm::TriggerResults>  handle_trigResults;
	iEvent.getByToken(token_trigResults, handle_trigResults);

	edm::TriggerNames triggerNames = iEvent.triggerNames(*handle_trigResults);

	for(unsigned int i_trig=0; i_trig<triggerNames.size(); ++i_trig) {
		std::string pathName = triggerNames.triggerName(i_trig);
		std::string str_isFired = "FAIL";
		if( handle_trigResults->accept(i_trig) ) str_isFired = "PASS";
		cout << pathName << ": " << str_isFired << endl;
	}

	// -- beam spot & vertex info.
	edm::Handle< reco::BeamSpot > beamSpot;
	iEvent.getByToken(beamSpot_, beamSpot);
	vertex_ = beamSpot->position();

	// // -- filter info.
	// edm::Handle<trigger::TriggerEvent> handle_trigEvent;
	// iEvent.getByToken(token_trigEvent, handle_trigEvent);

	// // -- hltIterL3MuonsNoID info.
	// edm::Handle< vector<reco::Muon> > handle_L3MuonsNoID;
	// if( iEvent.getByToken(token_L3MuonsNoID, handle_L3MuonsNoID) ) {
	// 	printMuonInfo( handle_L3MuonsNoID, "L3MuonsNoID" );
	// }
	// else
	// 	cout << "hltIterL3MuonsNoID: not reconstructed" << endl;

	// cout << endl;

	////////////////
	// -- Jets -- //
	////////////////
	edm::Handle< vector<reco::CaloJet> > caloJet;
	if( iEvent.getByToken(caloJet_, caloJet) )
	{
		cout << "========================================================" << endl;
		unsigned int nObj = caloJet->size();
		cout << "hltBSoftMuonDiJet40L1FastJetL25Jets (" << nObj << " objects)" << endl;

		for(unsigned int i_obj=0; i_obj<nObj; ++i_obj) {
			const reco::CaloJet& jet(caloJet->at(i_obj));
			printf("[%03d object] (pt, eta, phi) = (%9.3lf, %5.3lf, %5.3lf)\n", 
				       i_obj, jet.pt(), jet.eta(), jet.phi() );
		}

		cout << "========================================================" << endl;
		cout << endl;
	}
	else
		cout << "hltBSoftMuonDiJet40L1FastJetL25Jets: not reconstructed" << endl;

	//////////////////
	// -- tracks -- //
	//////////////////

	// -- Merged muon track info.
	edm::Handle< vector<reco::Track> > handle_muonMergedTrack;
	if( iEvent.getByToken(token_muonMergedTrack, handle_muonMergedTrack) ) {
		printRecoTrackInfo( handle_muonMergedTrack, "hltIterL3MuonAndMuonFromL1Merged");
	}
	else
		cout << "hltIterL3MuonAndMuonFromL1Merged: not reconstructed" << endl;

	// -- hltIterL3MuonTracks info.
	edm::Handle< vector<reco::Track> > handle_hltIterL3MuonTracks;
	if( iEvent.getByToken(token_hltIterL3MuonTracks, handle_hltIterL3MuonTracks) ) {
		printRecoTrackInfo( handle_hltIterL3MuonTracks, "hltIterL3MuonTracks");
	}
	else
		cout << "hltIterL3MuonTracks: not reconstructed" << endl;

	// -- hltBSoftMuonMu5L3 info
	edm::Handle< edm::View<reco::Track > > softMuTrackPt5;
	if( iEvent.getByToken(softMuTrackPt5_, softMuTrackPt5) )
	{
		cout << "========================================================" << endl;
		unsigned int nObj = softMuTrackPt5->size();
		cout << "hltBSoftMuonMu5L3 (" << nObj << " objects)" << endl;

		for(unsigned int i_obj=0; i_obj<nObj; ++i_obj) {
			const reco::Track& track(softMuTrackPt5->at(i_obj));
			printf("[%03d object] (pt, eta, phi, qualityMask, |dxy(BS)|, |dsz(BS)|, normChi2, nHit, nPixelHit, nLayer, n3DLayer, isGood) = \n(%9.3lf, %5.3lf, %5.3lf, %02d, %5.3lf, %5.3lf, %5.3lf, %02d, %02d, %02d, %02d, %02d)\n", 
				       i_obj, track.pt(), track.eta(), track.phi(), track.qualityMask(), fabs(track.dxy(vertex_)), fabs(track.dsz(vertex_)), track.normalizedChi2(),
				       track.hitPattern().numberOfValidHits(), track.hitPattern().numberOfValidPixelHits(), track.hitPattern().trackerLayersWithMeasurement(),
				       track.hitPattern().pixelLayersWithMeasurement() + track.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), isGoodTrack(track) );
		}

		cout << "========================================================" << endl;
		cout << endl;
	}
	else
		cout << "hltBSoftMuonMu5L3: not reconstructed" << endl;


	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
}

void MHPrintBTagMuTriggerInfo::printRecoTrackInfo(const edm::Handle< vector<reco::Track> >& handle_track, const std::string& type)
{
	cout << "========================================================" << endl;
	unsigned int nObj = handle_track->size();
	cout << type << " (" << nObj << " objects)" << endl;

	for(unsigned int i_obj=0; i_obj<nObj; ++i_obj)
	{
		const reco::Track& track(handle_track->at(i_obj));
		printf("[%03d object] (pt, eta, phi, qualityMask, |dxy(BS)|, |dsz(BS)|, normChi2, nHit, nPixelHit, nLayer, n3DLayer, isGood) = \n(%9.3lf, %5.3lf, %5.3lf, %02d, %5.3lf, %5.3lf, %5.3lf, %02d, %02d, %02d, %02d, %02d)\n", 
			       i_obj, track.pt(), track.eta(), track.phi(), track.qualityMask(), fabs(track.dxy(vertex_)), fabs(track.dsz(vertex_)), track.normalizedChi2(),
			       track.hitPattern().numberOfValidHits(), track.hitPattern().numberOfValidPixelHits(), track.hitPattern().trackerLayersWithMeasurement(),
			       track.hitPattern().pixelLayersWithMeasurement() + track.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), isGoodTrack(track) );
	}
	cout << "========================================================" << endl;
	cout << endl;
}

void MHPrintBTagMuTriggerInfo::printMuonInfo(const edm::Handle< vector<reco::Muon> >& handle_muon, const std::string& type)
{
	cout << "========================================================" << endl;
	unsigned int nObj = handle_muon->size();
	cout << type << " (" << nObj << " objects)" << endl;

	for(unsigned int i_mu=0; i_mu<nObj; ++i_mu) {
		const reco::Muon& muon(handle_muon->at(i_mu));
		printf("[%03d object] (pt, eta, phi) = (%9.3lf, %5.3lf, %5.3lf)\n", 
			       i_mu, muon.pt(), muon.eta(), muon.phi() );
	}

	cout << "========================================================" << endl;
	cout << endl;
}

bool MHPrintBTagMuTriggerInfo::isGoodTrack(const reco::Track& t)
{
	bool flag = false;

	double ptMin_ = 5.0;
	double minRapidity_ = -5.0;
	double maxRapidity_ = 5.0;
	double tip_ = 120.0;
	double lip_ = 300.0;
	double maxChi2_ = 10000.0;
	double minPhi = -3.2;
	double maxPhi = 3.2;
	int    minHit_ = 0;
	int    minPixelHit_ = 0;
	int    minLayer_ = 0;
	int    min3DLayer_ = 0;

	double meanPhi_ = (minPhi + maxPhi) / 2.0;
	double rangePhi_ = (maxPhi - minPhi) / 2.0;

	const auto dphi = deltaPhi(t.phi(), meanPhi_);

	if(  t.hitPattern().numberOfValidHits() >= minHit_ &&
       t.hitPattern().numberOfValidPixelHits() >= minPixelHit_ &&
       t.hitPattern().trackerLayersWithMeasurement() >= minLayer_ &&
       t.hitPattern().pixelLayersWithMeasurement() +
       t.hitPattern().numberOfValidStripLayersWithMonoAndStereo() >= min3DLayer_ &&
       fabs(t.pt()) >= ptMin_ &&
       t.eta() >= minRapidity_ && t.eta() <= maxRapidity_ &&
       dphi >= -rangePhi_ && dphi <= rangePhi_ &&
       fabs(t.dxy(vertex_)) <= tip_ &&
       fabs(t.dsz(vertex_)) <= lip_  &&
       t.normalizedChi2()<=maxChi2_ ) flag = true;

	return flag;
}

DEFINE_FWK_MODULE(MHPrintBTagMuTriggerInfo);