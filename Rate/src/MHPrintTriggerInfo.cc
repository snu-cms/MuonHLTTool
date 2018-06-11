// -- To investigate the increase of rates in MET triggers
// -- e.g. HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v8

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include <vector>
#include <iostream>

using namespace std;
using namespace edm;

class MHPrintTriggerInfo : public edm::EDAnalyzer {
public:
  MHPrintTriggerInfo(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

  void triggerFilterInfo (const edm::Handle < trigger::TriggerEvent >&, const std::string&);
  void printRecoTrackInfo(const edm::Handle < vector<reco::Track>   >&, const std::string&);
private:
	edm::EDGetTokenT< edm::TriggerResults          > token_trigResults;  
	edm::EDGetTokenT< trigger::TriggerEvent        > token_trigEvent;
	edm::EDGetTokenT< vector<reco::PFMET>          > token_PFMET;
	edm::EDGetTokenT< vector<reco::Track>          > token_muonMergedTrack;
	edm::EDGetTokenT< vector<reco::Muon>           > token_hltMuons;
	edm::EDGetTokenT< vector<reco::Muon>           > token_L3Muons;
	edm::EDGetTokenT< vector<reco::Track>          > token_hltPFMuonMerging;
	edm::EDGetTokenT< vector<reco::Track>          > token_hltMergedTracks;
	edm::EDGetTokenT< vector<reco::Muon>           > token_L3MuonsNoID;
	edm::EDGetTokenT< vector<float>                > token_MVAValue;
	edm::EDGetTokenT< vector<reco::Track>          > token_hltIterL3MuonTracks;
	// edm::EDGetTokenT< vector<reco::MuonTrackLinks> > token_hltMuonLinks;
	// edm::EDGetTokenT< vector<reco::MuonTrackLinks> > token_hltL3MuonsIterL3Links;
};


MHPrintTriggerInfo::MHPrintTriggerInfo(const edm::ParameterSet& iConfig):
token_trigResults           ( consumes< edm::TriggerResults          > ( iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")) ),
token_trigEvent             ( consumes< trigger::TriggerEvent        > ( iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"))   ),
token_PFMET                 ( consumes< vector<reco::PFMET>          > ( edm::InputTag("hltPFMETProducer::MYHLT"))                       ),
token_muonMergedTrack       ( consumes< vector<reco::Track>          > ( edm::InputTag("hltIterL3MuonAndMuonFromL1Merged::MYHLT"))       ),
token_hltMuons              ( consumes< vector<reco::Muon>           > ( edm::InputTag("hltMuons::MYHLT"))                               ),
token_L3Muons               ( consumes< vector<reco::Muon>           > ( edm::InputTag("hltIterL3Muons::MYHLT"))                         ),
token_hltPFMuonMerging      ( consumes< vector<reco::Track>          > ( edm::InputTag("hltPFMuonMerging::MYHLT"))                       ),
token_hltMergedTracks       ( consumes< vector<reco::Track>          > ( edm::InputTag("hltMergedTracks::MYHLT"))                        ),
token_L3MuonsNoID           ( consumes< vector<reco::Muon>           > ( edm::InputTag("hltIterL3MuonsNoID::MYHLT"))                     ),
token_MVAValue              ( consumes< vector<float>                > ( edm::InputTag("hltIterL3MuonAndMuonFromL1Merged", "MVAValues", "MYHLT")) ),
token_hltIterL3MuonTracks   ( consumes< vector<reco::Track>          > ( edm::InputTag("hltIterL3MuonTracks::MYHLT"))                    )
// token_hltMuonLinks          ( consumes< vector<reco::MuonTrackLinks> > ( edm::InputTag("hltMuonLinks::MYHLT"))                           ),
// token_hltL3MuonsIterL3Links ( consumes< vector<reco::MuonTrackLinks> > ( edm::InputTag("hltL3MuonsIterL3Links::MYHLT"))                  )
{

}

void MHPrintTriggerInfo::analyze( const edm::Event& iEvent, const edm::EventSetup& iEventSetup) {
	// -- print run:lumi:event number
	cout << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event() << endl;

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

	// -- filter info.
	edm::Handle<trigger::TriggerEvent> handle_trigEvent;
	iEvent.getByToken(token_trigEvent, handle_trigEvent);

	triggerFilterInfo(handle_trigEvent, "hltL3fL1DiMu3SQETM50f0PreFiltered3::MYHLT");
	triggerFilterInfo(handle_trigEvent, "hltDoubleMuon3Mass3p8to60Filter::MYHLT");
	triggerFilterInfo(handle_trigEvent, "hltDoubleMuon3Mass3p8to60DZFilter::MYHLT");
	triggerFilterInfo(handle_trigEvent, "hltMET40::MYHLT");
	triggerFilterInfo(handle_trigEvent, "hltPFMET90::MYHLT");
	triggerFilterInfo(handle_trigEvent, "hltPFMHTNoMuTightID90::MYHLT");

	// -- hltIterL3MuonsNoID info.
	edm::Handle< vector<reco::Muon> > handle_L3MuonsNoID;
	if( iEvent.getByToken(token_L3MuonsNoID, handle_L3MuonsNoID) ) {
		unsigned int nObj = handle_L3MuonsNoID->size();
		cout << "hltIterL3MuonsNoID (" << nObj << "objects)" << endl;

		for(unsigned int i_muon=0; i_muon<nObj; ++i_muon) {
			const reco::Muon& muon(handle_L3MuonsNoID->at(i_muon));
			cout << "\t(pt, eta, phi) = (" << muon.pt() << ", " << muon.eta() << ", " << muon.phi() << endl;
		}
	}
	else
		cout << "hltIterL3MuonsNoID: not reconstructed" << endl;

	cout << endl;

	// -- hltIterL3Muons info.
	edm::Handle< vector<reco::Muon> > handle_L3Muons;
	if( iEvent.getByToken(token_L3Muons, handle_L3Muons) ) {
		cout << "hltIterL3Muons" << endl;
		for(unsigned int i_muon=0; i_muon<handle_L3Muons->size(); ++i_muon) {
			const reco::Muon& muon(handle_L3Muons->at(i_muon));
			cout << "\t(pt, eta, phi) = (" << muon.pt() << ", " << muon.eta() << ", " << muon.phi() << ") => " << i_muon << "th object" << endl;
		}
	}
	else
		cout << "hltIterL3Muons: not reconstructed" << endl;

	cout << endl;

	// -- hltMergedTracks info.
	edm::Handle< vector<reco::Track> > handle_hltMergedTracks;
	if( iEvent.getByToken(token_hltMergedTracks, handle_hltMergedTracks) ) {
		printRecoTrackInfo( handle_hltMergedTracks, "hltMergedTracks");
	}
	else
		cout << "hltMergedTracks: not reconstructed" << endl;


	// -- Merged muon track info.
	edm::Handle< vector<reco::Track> > handle_muonMergedTrack;
	if( iEvent.getByToken(token_muonMergedTrack, handle_muonMergedTrack) ) {
		printRecoTrackInfo( handle_muonMergedTrack, "hltIterL3MuonAndMuonFromL1Merged");
	}
	else
		cout << "hltIterL3MuonAndMuonFromL1Merged: not reconstructed" << endl;

	cout << endl;


	// -- hltIterL3MuonTracks info.
	edm::Handle< vector<reco::Track> > handle_hltIterL3MuonTracks;
	if( iEvent.getByToken(token_hltIterL3MuonTracks, handle_hltIterL3MuonTracks) ) {
		printRecoTrackInfo( handle_hltIterL3MuonTracks, "hltIterL3MuonTracks");
	}
	else
		cout << "hltIterL3MuonTracks: not reconstructed" << endl;

	cout << endl;


	// -- hltPFMuonMerging info.
	edm::Handle< vector<reco::Track> > handle_hltPFMuonMerging;
	if( iEvent.getByToken(token_hltPFMuonMerging, handle_hltPFMuonMerging) ) {
		printRecoTrackInfo( handle_hltPFMuonMerging, "hltPFMuonMerging");
	}
	else
		cout << "hltPFMuonMerging: not reconstructed" << endl;

	cout << endl;

	// -- hltMuons info.
	edm::Handle< vector<reco::Muon> > handle_hltMuons;
	if( iEvent.getByToken(token_hltMuons, handle_hltMuons) ) {
		cout << "hltMuons" << endl;
		for(unsigned int i_muon=0; i_muon<handle_hltMuons->size(); ++i_muon) {
			const reco::Muon& muon(handle_hltMuons->at(i_muon));
			cout << "\t(pt, eta, phi) = (" << muon.pt() << ", " << muon.eta() << ", " << muon.phi() << ") => " << i_muon << "th object" << endl;
		}
	}
	else
		cout << "hltMuons: not reconstructed" << endl;

	cout << endl;

	// -- PFMET info.
	edm::Handle<vector<reco::PFMET>> handle_PFMET;
	if( iEvent.getByToken(token_PFMET, handle_PFMET) ) {
		const reco::PFMET& pfMET(handle_PFMET->front());
		cout << "PFMET: " << pfMET.pt() << endl;
		cout << "\tSumEt: " << pfMET.sumEt() << endl;
		cout << "\tMuonEt: " << pfMET.MuonEt() << endl;
		cout << "\tNeutralEMEt: " << pfMET.NeutralEMEt() << endl;
		cout << "\tNeutralHadEt: " << pfMET.NeutralHadEt() << endl;
		cout << "\tChargedEMEt: " << pfMET.ChargedEMEt() << endl;
		cout << "\tChargedHadEt: " << pfMET.ChargedHadEt() << endl;
		cout << "\tType6Et: " << pfMET.Type6Et() << endl;
		cout << "\tType7Et: " << pfMET.Type7Et() << endl;
		cout << "\tType7Et: " << pfMET.Type7Et() << endl;
	}
	else
		cout << "PFMET is not reconstructed" << endl;

	cout << "\nMVA value of hltIterL3MuonAndMuonFromL1Merged" << endl;
	edm::Handle< vector<float> > handle_MVAValue;
	if( iEvent.getByToken(token_MVAValue, handle_MVAValue) ) {
		for(unsigned int i=0; i<(*handle_MVAValue).size(); i++) {
			cout << i << "th MVA value: " << (*handle_MVAValue)[i] << endl;
		}
	}
	else
		cout << "Fail to load MVA value" << endl;


	cout << endl;
	cout << endl;
}

void MHPrintTriggerInfo::triggerFilterInfo(const edm::Handle<trigger::TriggerEvent>& handle_trigEvent, const std::string& filterName) {
  
  const trigger::size_type nFilters(handle_trigEvent->sizeFilters());
  for( trigger::size_type i_filter=0; i_filter<nFilters; i_filter++) {
  	std::string filterTag = handle_trigEvent->filterTag(i_filter).encode();

  	if( filterTag == filterName ) {
  		trigger::Keys objectKeys = handle_trigEvent->filterKeys(i_filter);
  		Int_t nObj = objectKeys.size();
  		cout << filterName <<  ": # objects = " << nObj	<< endl;

  		const trigger::TriggerObjectCollection& triggerObjects(handle_trigEvent->getObjects());

  		for( trigger::size_type i_key=0; i_key<nObj; i_key++) {
  		  trigger::size_type objKey = objectKeys.at(i_key);
  		  const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
  		  cout << "\t(pt, eta, phi) = (" << triggerObj.pt() << ", " << triggerObj.eta() << ", " << triggerObj.phi() << ") => " << i_key << "th object" << endl;
  		}
  	}

  }
}

void MHPrintTriggerInfo::printRecoTrackInfo(const edm::Handle< vector<reco::Track> >& handle_track, const std::string& type) {
	unsigned int nObj = handle_track->size();
	cout << type << "(# object = " << nObj << ")" << endl;

	for(unsigned int i_track=0; i_track<nObj; ++i_track) {
		const reco::Track& track(handle_track->at(i_track));
		cout << "\t(pt, eta, phi, qualityMask) = (" \
		     << track.pt()          << ", " \
		     << track.eta()         << ", " \
		     << track.phi()         << ", " \
		     << track.qualityMask() << ") " << endl;
	}
}

DEFINE_FWK_MODULE(MHPrintTriggerInfo);