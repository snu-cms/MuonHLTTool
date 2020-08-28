#ifndef MuonHLTobjCorrelator_h
#define MuonHLTobjCorrelator_h 1

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

namespace MuonHLTobjCorrelator {
  class L1TTTrack {
  private:
    float pt;
    float eta;
    float phi;
    float z0;
    float rInv;

  public:
    L1TTTrack(const TTTrack< Ref_Phase2TrackerDigi_ > theTrack) {
      pt = static_cast<float>( theTrack.momentum().perp() );
      eta = static_cast<float>( theTrack.momentum().eta() );
      phi = static_cast<float>( theTrack.momentum().phi() );
      z0 = static_cast<float>( theTrack.z0() );
      rInv = static_cast<float>( theTrack.rInv() );
    }

    bool operator<(const L1TTTrack& other) const {
      bool checkSequence = (this->pt!=other.pt);
      if (checkSequence) return (this->pt < other.pt);

      checkSequence = (this->eta!=other.eta);
      if (checkSequence) return (this->eta < other.eta);

      checkSequence = (this->phi!=other.phi);
      if (checkSequence) return (this->phi < other.phi);

      checkSequence = (this->z0!=other.z0);
      if (checkSequence) return (this->z0 < other.z0);

      return (this->rInv < other.rInv);
    }
  };
}

#endif
