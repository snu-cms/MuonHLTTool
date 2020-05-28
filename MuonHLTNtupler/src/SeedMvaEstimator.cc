#include "MuonHLTTool/MuonHLTNtupler/interface/SeedMvaEstimator.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"

using namespace std;

SeedMvaEstimator::SeedMvaEstimator(const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std) {
  gbrForest_  = createGBRForest(weightsfile);
  scale_mean_ = scale_mean;
  scale_std_  = scale_std;
}

SeedMvaEstimator::~SeedMvaEstimator() {}

namespace {
  enum inputIndexes {
    kTsosErr0,         // 0
    kTsosErr1,         // 1
    kTsosErr2,         // 2
    kTsosErr3,         // 3
    kTsosErr4,         // 4
    kTsosErr5,         // 5
    kTsosErr6,         // 6
    kTsosErr7,         // 7
    kTsosErr8,         // 8
    kTsosErr9,         // 9
    kTsosErr10,        // 10
    kTsosErr11,        // 11
    kTsosErr12,        // 12
    kTsosErr13,        // 13
    kTsosErr14,        // 14
    kTsosDxdz,         // 15
    kTsosDydz,         // 16
    kTsosQbp,          // 17
    kTsosCharge,       // 18
    kDRdRL1SeedP,      // 19
    kDPhidRL1SeedP,    // 20
    kDRdPhiL1SeedX,    // 21
    kDPhidPhiL1SeedX,  // 22
    kDRdRL2SeedP,      // 23
    kDPhidRL2SeedP,    // 24
    kDRdPhiL2SeedX,    // 25
    kDPhidPhiL2SeedX,  // 26
    kDRL1TkMu,         // 27
    kDPhiL1TkMu,       // 28
    kLast              // 29
  };
}  // namespace

void SeedMvaEstimator::getL1MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  float& dRdRL1SeedP,
  float& dPhidRL1SeedP,
  float& dRdPhiL1SeedX,
  float& dPhidPhiL1SeedX ) const {

  for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
  {
    if(ibx != 0) continue; // -- only take when ibx == 0 -- //
    for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
    {
      l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

      // FIXME: 7 should be configurable
      if(ref_L1Mu->hwQual() < 7)
        continue;

      float dR_L1SeedP_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_p.eta(), global_p.phi());
      float dPhi_L1SeedP_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_p.phi());
      float dR_L1SeedX_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_x.eta(), global_x.phi());
      float dPhi_L1SeedX_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_x.phi());

      if( dR_L1SeedP_AtVtx < dRdRL1SeedP ) {
        dRdRL1SeedP = dR_L1SeedP_AtVtx;
        dPhidRL1SeedP = dPhi_L1SeedP_AtVtx;
      }
      if( fabs(dPhi_L1SeedX_AtVtx) < fabs(dPhidPhiL1SeedX) ) {
        dRdPhiL1SeedX = dR_L1SeedX_AtVtx;
        dPhidPhiL1SeedX = dPhi_L1SeedX_AtVtx;
      }
    }
  }
}

void SeedMvaEstimator::getL2MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  float& dRdRL2SeedP,
  float& dPhidRL2SeedP,
  float& dRdPhiL2SeedX,
  float& dPhidPhiL2SeedX ) const {

  for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
  {
    reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, i_L2);

    float dR_L2SeedP   = reco::deltaR( *ref_L2Mu, global_p);
    float dPhi_L2SeedP = reco::deltaPhi( ref_L2Mu->phi(), global_p.phi());
    float dR_L2SeedX   = reco::deltaR( *ref_L2Mu, global_x);
    float dPhi_L2SeedX = reco::deltaPhi( ref_L2Mu->phi(), global_x.phi());

    if( dR_L2SeedP < dRdRL2SeedP ) {
      dRdRL2SeedP = dR_L2SeedP;
      dPhidRL2SeedP = dPhi_L2SeedP;
    }
    if( fabs(dPhi_L2SeedX) < fabs(dPhidPhiL2SeedX) ) {
      dRdPhiL2SeedX = dR_L2SeedX;
      dPhidPhiL2SeedX = dPhi_L2SeedX;
    }
  }
}

void SeedMvaEstimator::getL1TTVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu,
  float& DRL1TkMu,
  float& DPhiL1TkMu ) const {

  for(auto L1TkMu=h_L1TkMu->begin(); L1TkMu!=h_L1TkMu->end(); ++L1TkMu)
  {
    auto TkRef = L1TkMu->trkPtr();
    float DRL1TkMu_tmp   = reco::deltaR( *TkRef, global_p);
    float DPhiL1TkMu_tmp = reco::deltaPhi( TkRef->phi(), global_p.phi());
    if( DRL1TkMu_tmp < DRL1TkMu ) {
      DRL1TkMu   = DRL1TkMu_tmp;
      DPhiL1TkMu = DPhiL1TkMu_tmp;
    }
  }
}

float SeedMvaEstimator::computeMva( const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu
) const {

  float var[kLast]{};

  var[kTsosErr0]   = seed.startingState().error(0);
  var[kTsosErr1]   = seed.startingState().error(1);
  var[kTsosErr2]   = seed.startingState().error(2);
  var[kTsosErr3]   = seed.startingState().error(3);
  var[kTsosErr4]   = seed.startingState().error(4);
  var[kTsosErr5]   = seed.startingState().error(5);
  var[kTsosErr6]   = seed.startingState().error(6);
  var[kTsosErr7]   = seed.startingState().error(7);
  var[kTsosErr8]   = seed.startingState().error(8);
  var[kTsosErr9]   = seed.startingState().error(9);
  var[kTsosErr10]  = seed.startingState().error(10);
  var[kTsosErr11]  = seed.startingState().error(11);
  var[kTsosErr12]  = seed.startingState().error(12);
  var[kTsosErr13]  = seed.startingState().error(13);
  var[kTsosErr14]  = seed.startingState().error(14);
  var[kTsosDxdz]   = seed.startingState().parameters().dxdz();
  var[kTsosDydz]   = seed.startingState().parameters().dydz();
  var[kTsosQbp]    = seed.startingState().parameters().qbp();
  var[kTsosCharge] = seed.startingState().parameters().charge();

  // FIXME: should be configurable
  float initDRdPhi = 99999.;

  float dRdRL1SeedP = initDRdPhi;
  float dPhidRL1SeedP = initDRdPhi;
  float dRdPhiL1SeedX = initDRdPhi;
  float dPhidPhiL1SeedX = initDRdPhi;
  getL1MuonVariables( seed, global_p, global_x, h_L1Muon, dRdRL1SeedP, dPhidRL1SeedP, dRdPhiL1SeedX, dPhidPhiL1SeedX );

  float dRdRL2SeedP = initDRdPhi;
  float dPhidRL2SeedP = initDRdPhi;
  float dRdPhiL2SeedX = initDRdPhi;
  float dPhidPhiL2SeedX = initDRdPhi;
  getL2MuonVariables( seed, global_p, global_x, h_L2Muon, dRdRL2SeedP, dPhidRL2SeedP, dRdPhiL2SeedX, dPhidPhiL2SeedX );

  float DRL1TkMu = initDRdPhi;
  float DPhiL1TkMu = initDRdPhi;
  getL1TTVariables( seed, global_p, global_x, h_L1TkMu, DRL1TkMu, DPhiL1TkMu );

  var[kDRdRL1SeedP]     = dRdRL1SeedP;
  var[kDPhidRL1SeedP]   = dPhidRL1SeedP;
  var[kDRdPhiL1SeedX]   = dRdPhiL1SeedX;
  var[kDPhidPhiL1SeedX] = dPhidPhiL1SeedX;
  var[kDRdRL2SeedP]     = dRdRL2SeedP;
  var[kDPhidRL2SeedP]   = dPhidRL2SeedP;
  var[kDRdPhiL2SeedX]   = dRdPhiL2SeedX;
  var[kDPhidPhiL2SeedX] = dPhidPhiL2SeedX;
  var[kDRL1TkMu]        = DRL1TkMu;
  var[kDPhiL1TkMu]      = DPhiL1TkMu;

  for(int iv=0; iv<kLast; ++iv) {
    var[iv] = (var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
  }

  return gbrForest_->GetResponse( var );
}
