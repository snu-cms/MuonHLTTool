#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <vector>
#include <iostream>

using namespace std;
using namespace edm;

class MHFindTriggerFiredEvents : public edm::EDAnalyzer {
  public:
  	explicit MHFindTriggerFiredEvents(const edm::ParameterSet&);
  	void analyze(const edm::Event&, const edm::EventSetup&);
  private:
  	edm::EDGetTokenT<edm::TriggerResults> token_trigResults;
  	string trigName;
};

MHFindTriggerFiredEvents::MHFindTriggerFiredEvents(const edm::ParameterSet& iConfig):
token_trigResults( consumes< edm::TriggerResults > (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults")) ),
trigName( iConfig.getUntrackedParameter<std::string>("trigName") )
{

}

void MHFindTriggerFiredEvents::analyze( const edm::Event& iEvent, const edm::EventSetup& iEventSetup) {
  edm::Handle<edm::TriggerResults>  handle_trigResults;
  iEvent.getByToken(token_trigResults, handle_trigResults);

  edm::TriggerNames triggerNames = iEvent.triggerNames(*handle_trigResults);

  for(unsigned int i_trig=0; i_trig<triggerNames.size(); ++i_trig) {

    // -- if accepted
    if( handle_trigResults->accept(i_trig) ) {
      std::string pathName = triggerNames.triggerName(i_trig);

      // -- if the trigger is fired in this event
      if( pathName.find( this->trigName ) != std::string::npos ) {
        // cout << "pathName: " << pathName << " is found!" << endl;
        cout << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event() << endl;
      }
    }

  }
}

DEFINE_FWK_MODULE(MHFindTriggerFiredEvents);