#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


class DYmuSkimmer : public edm::EDFilter {
public:
  explicit DYmuSkimmer(const edm::ParameterSet&);

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >  t_pu_info;
  edm::EDGetTokenT< reco::GenParticleCollection >     t_gen_particle;
  double        _pu_min;
  double        _pu_max;
};


DYmuSkimmer::DYmuSkimmer(const edm::ParameterSet& iConfig)
  : t_pu_info( consumes< std::vector<PileupSummaryInfo> >(   iConfig.getParameter<edm::InputTag>("pu_info_src") ) ),
    t_gen_particle( consumes< reco::GenParticleCollection >( iConfig.getParameter<edm::InputTag>("gen_particle_src") ) ),
    _pu_min(      iConfig.getParameter< double >("pu_min")),
    _pu_max(      iConfig.getParameter< double >("pu_max"))
{
}

bool DYmuSkimmer::filter(edm::Event& iEvent, const edm::EventSetup&) {

  // -- PU filter -- //
  bool pass_pu = false;
  double pu = -1;
  edm::Handle<std::vector< PileupSummaryInfo > > h_pu_info;
  if( iEvent.getByToken(t_pu_info,h_pu_info) ) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = h_pu_info->begin(); PVI != h_pu_info->end(); ++PVI){
      if(PVI->getBunchCrossing()==0)
      {
        pu = PVI->getTrueNumInteractions();
        continue;
      }
    }
  }

  if( pu > _pu_min && pu < _pu_max )
    pass_pu = true;

  // -- Gen mumu filter -- //
  bool found0 = false;
  bool found1 = false;
  edm::Handle<reco::GenParticleCollection> h_gen_particle;
  if( iEvent.getByToken(t_gen_particle,h_gen_particle) ) {
    reco::GenParticleCollection::const_iterator genp = h_gen_particle->begin();
    for (; genp != h_gen_particle->end(); genp++) {
      if(genp->pdgId() ==  13 && genp->isHardProcess()==1) {
        // const reco::Candidate* m = genp->mother();
        // if(m->pdgId()==23){
          found0=true;
        // }
      }

      if(genp->pdgId()==-13 && genp->isHardProcess()==1) {
        // const reco::Candidate* m2 = genp->mother();
        // if(m2->pdgId()==23){
          found1=true;
        // }
      }
    }
  }

  return (pass_pu && found0 && found1); 
}

DEFINE_FWK_MODULE(DYmuSkimmer);
