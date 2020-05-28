#ifndef __PhysicsTools_PatAlgos_SeedMvaEstimator__
#define __PhysicsTools_PatAlgos_SeedMvaEstimator__

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"

#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include <memory>
#include <string>

class GBRForest;

namespace edm {
  class FileInPath;
}

class SeedMvaEstimator {
public:
  SeedMvaEstimator( const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std );
  ~SeedMvaEstimator();

  float computeMva( const TrajectorySeed&,
    GlobalVector,
    GlobalPoint,
    edm::Handle<l1t::MuonBxCollection>,
    edm::Handle<reco::RecoChargedCandidateCollection>,
    edm::Handle<l1t::TkMuonCollection>
  ) const;

private:
  std::unique_ptr<const GBRForest> gbrForest_;

  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;

  void getL1MuonVariables( const TrajectorySeed&, GlobalVector, GlobalPoint, edm::Handle<l1t::MuonBxCollection>, float&, float&, float&, float& ) const;
  void getL2MuonVariables( const TrajectorySeed&, GlobalVector, GlobalPoint, edm::Handle<reco::RecoChargedCandidateCollection>, float&, float&, float&, float& ) const;
  void getL1TTVariables(   const TrajectorySeed&, GlobalVector, GlobalPoint, edm::Handle<l1t::TkMuonCollection>, float&, float& ) const;
};

#endif
