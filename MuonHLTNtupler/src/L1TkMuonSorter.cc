// -*- C++ -*-
//
// Package:    tmp/L1TkMuonSorter
// Class:      L1TkMuonSorter
// 
/**\class L1TkMuonSorter L1TkMuonSorter.cc tmp/L1TkMuonSorter/plugins/L1TkMuonSorter.cc

 Description: [one line class summary]

 Implementation:
	  [Notes on implementation]
*/
//
// Original Author:  Minseok Oh
//         Created:  Wed, 13 May 2020 08:42:57 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

//
// class declaration
//
class L1TkMuonSorter : public edm::stream::EDProducer<> {
	public:
		explicit L1TkMuonSorter(const edm::ParameterSet&);
		~L1TkMuonSorter();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginStream(edm::StreamID) override;
		virtual void produce(edm::Event&, const edm::EventSetup&) override;
		virtual void endStream() override;

		// ----------member data ---------------------------
		const edm::EDGetTokenT<l1t::TkMuonCollection> t_L1TkMuon_;

		struct pt_sort {
			bool operator()(const l1t::TkMuon& lhs, const l1t::TkMuon& rhs) {
				return lhs.pt() > rhs.pt();
			}
		};
};

//
// constructors and destructor
//
L1TkMuonSorter::L1TkMuonSorter(const edm::ParameterSet& iConfig):
	t_L1TkMuon_(consumes<l1t::TkMuonCollection>(iConfig.getParameter<edm::InputTag>("L1TkMuons")))
{
	produces<l1t::TkMuonCollection>();
}

L1TkMuonSorter::~L1TkMuonSorter(){}


//
// member functions
//

// ------------ method called to produce the data  ------------
void L1TkMuonSorter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	auto output = std::make_unique<l1t::TkMuonCollection>();

	edm::Handle<l1t::TkMuonCollection> h_L1TkMu;
	iEvent.getByToken(t_L1TkMuon_, h_L1TkMu);

	for (unsigned int i = 0; i < h_L1TkMu->size(); ++i) {
		const l1t::TkMuon& L1TkMu(h_L1TkMu->at(i));
		output->push_back(L1TkMu);
	}

	sort(output->begin(), output->end(), pt_sort());

	iEvent.put(std::move(output));
}




// ------------ method called once each stream before processing any runs, lumis or events  ------------
void L1TkMuonSorter::beginStream(edm::StreamID){}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void L1TkMuonSorter::endStream(){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkMuonSorter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonSorter);
