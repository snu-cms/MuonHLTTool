#include "TH1D.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <vector>

using namespace std;
using namespace edm;

class MHMETAnalyzer : public edm::EDAnalyzer {
  public:
  	explicit MHMETAnalyzer(const edm::ParameterSet&);
  	void analyze(const edm::Event&, const edm::EventSetup&);
  private:
  	edm::EDGetTokenT<reco::PFMETCollection> token_onlinePFMET;
  	edm::EDGetTokenT<vector<reco::GenMET>> token_genMET;

  	TH1D* h_genMET;
  	TH1D* h_onlinePFMET;
  	TH1D* h_relDiffMET;
};

MHMETAnalyzer::MHMETAnalyzer(const edm::ParameterSet& iConfig):
token_onlinePFMET( consumes< reco::PFMETCollection > (iConfig.getUntrackedParameter<edm::InputTag>("onlinePFMET")) ),
token_genMET     ( consumes< vector<reco::GenMET>  > (iConfig.getUntrackedParameter<edm::InputTag>("genMET"     )) )
{
	edm::Service<TFileService> fs;
	this->h_genMET = fs->make<TH1D>("h_genMET", "", 10000, 0, 10000);
	this->h_onlinePFMET = fs->make<TH1D>("h_onlinePFMET", "", 10000, 0, 10000);
	this->h_relDiffMET = fs->make<TH1D>("h_relDiffMET", "", 2000, -10, 10);
}

void MHMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup) {
	// -- genMET
	double genMET = 9999;

	edm::Handle< vector<reco::GenMET> > Handle_genMET;
	if( iEvent.getByToken(token_genMET, Handle_genMET) )
		genMET = Handle_genMET->front().pt();

	// cout << "genMET: " << genMET << endl;
	h_genMET->Fill( genMET );

	// -- online PF MET
	double onlinePFMET = 9999;
	edm::Handle< reco::PFMETCollection > handle_onlinePFMET;
	if( iEvent.getByToken( token_onlinePFMET, handle_onlinePFMET ) )
		onlinePFMET = handle_onlinePFMET->front().pt();

	// cout << "onlinePFMET: " << onlinePFMET << endl;
	h_onlinePFMET->Fill( onlinePFMET );

	// -- relative difference between genMET and onlinePFMET
	double relDiff = 0;
	if( genMET != 0 )
		relDiff = (onlinePFMET - genMET ) / genMET;
	else
		relDiff = 9.9;

	h_relDiffMET->Fill( relDiff );
}

DEFINE_FWK_MODULE(MHMETAnalyzer);