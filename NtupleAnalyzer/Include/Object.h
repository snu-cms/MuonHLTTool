#pragma once
#include <TLorentzVector.h>

//customized header files
#include <Include/NtupleHandle.h>

namespace MuonHLT
{

const Double_t M_mu   = 0.1056583715;  // -- GeV -- //
const Double_t M_elec = 0.000510998; // -- GeV -- //
const Double_t M_tau  = 1.77682;      // -- GeV -- //

class Object
{
public:
  Double_t pt;
  Double_t eta;
  Double_t phi;
  Double_t mass;
  TLorentzVector vecP;

  Object() {}
};

class GenParticle : public Object
{
public:
  Int_t ID;
  Int_t status;
  Int_t mother;

  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t energy;
  Double_t charge;

  Int_t isPrompt;
  Int_t isPromptFinalState;
  Int_t isTauDecayProduct;
  Int_t isPromptTauDecayProduct;
  Int_t isDirectPromptTauDecayProductFinalState;
  Int_t isHardProcess;
  Int_t isLastCopy;
  Int_t isLastCopyBeforeFSR;
  Int_t isPromptDecayed;
  Int_t isDecayedLeptonHadron;
  Int_t fromHardProcessBeforeFSR;
  Int_t fromHardProcessDecayed;
  Int_t fromHardProcessFinalState;
  Int_t isMostlyLikePythia6Status3;

  GenParticle()
  {
    Init();
  }

  GenParticle(NtupleHandle* ntuple, Int_t index): GenParticle()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    ID     = ntuple->genParticle_ID[index];
    status = ntuple->genParticle_status[index];
    mother = ntuple->genParticle_mother[index];

    pt     = ntuple->genParticle_pt[index];
    eta    = ntuple->genParticle_eta[index];
    phi    = ntuple->genParticle_phi[index];
    px     = ntuple->genParticle_px[index];
    py     = ntuple->genParticle_py[index];
    pz     = ntuple->genParticle_pz[index];
    energy = ntuple->genParticle_energy[index];
    charge = ntuple->genParticle_charge[index];

    vecP.SetPxPyPzE( px, py, pz, energy );
    mass = vecP.M();

    isPrompt                                = ntuple->genParticle_isPrompt[index];
    isPromptFinalState                      = ntuple->genParticle_isPromptFinalState[index];
    isTauDecayProduct                       = ntuple->genParticle_isTauDecayProduct[index];
    isPromptTauDecayProduct                 = ntuple->genParticle_isPromptTauDecayProduct[index];
    isDirectPromptTauDecayProductFinalState = ntuple->genParticle_isDirectPromptTauDecayProductFinalState[index];
    isHardProcess                            = ntuple->genParticle_isHardProcess[index];
    isLastCopy                               = ntuple->genParticle_isLastCopy[index];
    isLastCopyBeforeFSR                      = ntuple->genParticle_isLastCopyBeforeFSR[index];
    isPromptDecayed                          = ntuple->genParticle_isPromptDecayed[index];
    isDecayedLeptonHadron                     = ntuple->genParticle_isDecayedLeptonHadron[index];
    fromHardProcessBeforeFSR                 = ntuple->genParticle_fromHardProcessBeforeFSR[index];
    fromHardProcessDecayed                   = ntuple->genParticle_fromHardProcessDecayed[index];
    fromHardProcessFinalState                = ntuple->genParticle_fromHardProcessFinalState[index];
    isMostlyLikePythia6Status3               = ntuple->genParticle_isMostlyLikePythia6Status3[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;

    ID = -999;
    status = -999;
    mother = -999;

    px = -999;
    py = -999;
    pz = -999;
    energy = -999;
    charge = -999;

    isPrompt = -999;
    isPromptFinalState = -999;
    isTauDecayProduct = -999;
    isPromptTauDecayProduct = -999;
    isDirectPromptTauDecayProductFinalState = -999;
    isHardProcess = -999;
    isLastCopy = -999;
    isLastCopyBeforeFSR = -999;
    isPromptDecayed = -999;
    isDecayedLeptonHadron = -999;
    fromHardProcessBeforeFSR = -999;
    fromHardProcessDecayed = -999;
    fromHardProcessFinalState = -999;
    isMostlyLikePythia6Status3 = -999;
  }
};

class Muon : public Object
{
public:
  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t energy;
  Double_t dB;
  Double_t charge;
  Int_t    isGLB;
  Int_t    isSTA;
  Int_t    isTRK;
  Int_t    isPF;
  Int_t    isTight;
  Int_t    isMedium;
  Int_t    isLoose;
  Int_t    isHighPt;
  Int_t    isSoft;
  Double_t iso03_sumPt;
  Double_t iso03_hadEt;
  Double_t iso03_emEt;
  Double_t PFIso03_charged;
  Double_t PFIso03_neutral;
  Double_t PFIso03_photon;
  Double_t PFIso03_sumPU;
  Double_t PFIso04_charged;
  Double_t PFIso04_neutral;
  Double_t PFIso04_photon;
  Double_t PFIso04_sumPU;
  Double_t PFCluster03_ECAL;
  Double_t PFCluster03_HCAL;
  Double_t PFCluster04_ECAL;
  Double_t PFCluster04_HCAL;
  Double_t normChi2_global;
  Int_t    nTrackerHit_global;
  Int_t    nTrackerLayer_global;
  Int_t    nPixelHit_global;
  Int_t    nMuonHit_global;
  Double_t normChi2_inner;
  Int_t    nTrackerHit_inner;
  Int_t    nTrackerLayer_inner;
  Int_t    nPixelHit_inner;
  Double_t pt_tuneP;
  Double_t ptError_tuneP;
  Double_t dxyVTX_best;
  Double_t dzVTX_best;
  Int_t    nMatchedStation;
  Int_t    nMatchedRPCLayer;
  Int_t    stationMask;

  Double_t relPFIso_dBeta;
  Double_t relTrkIso;

  Muon()
  {
    Init();
  }

  Muon(NtupleHandle* ntuple, Int_t index): Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt     = ntuple->muon_pt[index];
    eta    = ntuple->muon_eta[index];
    phi    = ntuple->muon_phi[index];
    charge = ntuple->muon_charge[index];

    px  = ntuple->muon_px[index];
    py  = ntuple->muon_py[index];
    pz  = ntuple->muon_pz[index];
    mass = MuonHLT::M_mu;
    energy = sqrt(px*px + py*py + pz*pz + mass*mass);
    vecP.SetPxPyPzE(px, py, pz, energy);

    isGLB = ntuple->muon_isGLB[index];
    isSTA = ntuple->muon_isSTA[index];
    isTRK = ntuple->muon_isTRK[index];
    isPF  = ntuple->muon_isPF[index];

    isTight  = ntuple->muon_isTight[index];
    isMedium = ntuple->muon_isMedium[index];
    isLoose  = ntuple->muon_isLoose[index];
    isHighPt = ntuple->muon_isHighPt[index];
    isSoft   = ntuple->muon_isSoft[index];

    iso03_sumPt = ntuple->muon_iso03_sumPt[index];
    iso03_hadEt = ntuple->muon_iso03_hadEt[index];
    iso03_emEt  = ntuple->muon_iso03_emEt[index];

    PFIso03_charged = ntuple->muon_PFIso03_charged[index];
    PFIso03_neutral = ntuple->muon_PFIso03_neutral[index];
    PFIso03_photon  = ntuple->muon_PFIso03_photon[index];
    PFIso03_sumPU   = ntuple->muon_PFIso03_sumPU[index];
    PFIso04_charged = ntuple->muon_PFIso04_charged[index];
    PFIso04_neutral = ntuple->muon_PFIso04_neutral[index];
    PFIso04_photon  = ntuple->muon_PFIso04_photon[index];
    PFIso04_sumPU   = ntuple->muon_PFIso04_sumPU[index];

    PFCluster03_ECAL = ntuple->muon_PFCluster03_ECAL[index];
    PFCluster03_HCAL = ntuple->muon_PFCluster03_HCAL[index];
    PFCluster04_ECAL = ntuple->muon_PFCluster04_ECAL[index];
    PFCluster04_HCAL = ntuple->muon_PFCluster04_HCAL[index];

    normChi2_global      = ntuple->muon_normChi2_global[index];
    nTrackerHit_global   = ntuple->muon_nTrackerHit_global[index];
    nTrackerLayer_global = ntuple->muon_nTrackerLayer_global[index];
    nPixelHit_global     = ntuple->muon_nPixelHit_global[index];
    nMuonHit_global      = ntuple->muon_nMuonHit_global[index];

    normChi2_inner      = ntuple->muon_normChi2_inner[index];
    nTrackerHit_inner   = ntuple->muon_nTrackerHit_inner[index];
    nTrackerLayer_inner = ntuple->muon_nTrackerLayer_inner[index];
    nPixelHit_inner     = ntuple->muon_nPixelHit_inner[index];

    pt_tuneP      = ntuple->muon_pt_tuneP[index];
    ptError_tuneP = ntuple->muon_ptError_tuneP[index];

    dxyVTX_best = ntuple->muon_dxyVTX_best[index];
    dzVTX_best  = ntuple->muon_dzVTX_best[index];
    dB          = ntuple->muon_dB[index];

    nMatchedStation  = ntuple->muon_nMatchedStation[index];
    nMatchedRPCLayer = ntuple->muon_nMatchedRPCLayer[index];
    stationMask      = ntuple->muon_stationMask[index];

    relPFIso_dBeta = (PFIso04_charged + max(0., PFIso04_neutral + PFIso04_photon - 0.5*PFIso04_sumPU))/pt;
    relTrkIso = iso03_sumPt / pt;
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    px = -999;
    py = -999;
    pz = -999;
    energy = -999;
    dB = -999;
    charge = -999;
    isGLB = -999;
    isSTA = -999;
    isTRK = -999;
    isPF = -999;
    isTight = -999;
    isMedium = -999;
    isLoose = -999;
    isHighPt = -999;
    isSoft = -999;
    iso03_sumPt = -999;
    iso03_hadEt = -999;
    iso03_emEt = -999;
    PFIso03_charged = -999;
    PFIso03_neutral = -999;
    PFIso03_photon = -999;
    PFIso03_sumPU = -999;
    PFIso04_charged = -999;
    PFIso04_neutral = -999;
    PFIso04_photon = -999;
    PFIso04_sumPU = -999;
    PFCluster03_ECAL = -999;
    PFCluster03_HCAL = -999;
    PFCluster04_ECAL = -999;
    PFCluster04_HCAL = -999;
    normChi2_global = -999;
    nTrackerHit_global = -999;
    nTrackerLayer_global = -999;
    nPixelHit_global = -999;
    nMuonHit_global = -999;
    normChi2_inner = -999;
    nTrackerHit_inner = -999;
    nTrackerLayer_inner = -999;
    nPixelHit_inner = -999;
    pt_tuneP = -999;
    ptError_tuneP = -999;
    dxyVTX_best = -999;
    dzVTX_best = -999;
    nMatchedStation = -999;
    nMatchedRPCLayer = -999;
    stationMask = -999;

    relPFIso_dBeta = -999;
    relTrkIso = -999;
  }
};

class L3Muon : public Object
{
public:
  Double_t charge;
  Double_t trkPt;

  L3Muon()
  {
    Init();
  }

  L3Muon(NtupleHandle* ntuple, Int_t index): L3Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt   = ntuple->L3Muon_pt[index];
    eta  = ntuple->L3Muon_eta[index];
    phi  = ntuple->L3Muon_phi[index];
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    charge = ntuple->L3Muon_charge[index];
    trkPt  = ntuple->L3Muon_trkPt[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    charge = -999;
    trkPt = -999;
  }
};

class L2Muon : public Object
{
public:
  Double_t charge;
  Double_t trkPt;

  L2Muon()
  {
    Init();
  }

  L2Muon(NtupleHandle* ntuple, Int_t index): L2Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt   = ntuple->L2Muon_pt[index];
    eta  = ntuple->L2Muon_eta[index];
    phi  = ntuple->L2Muon_phi[index];
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    charge = ntuple->L2Muon_charge[index];
    trkPt  = ntuple->L2Muon_trkPt[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    charge = -999;
    trkPt = -999;
  }
};

class TkMuon : public Object
{
public:
  Double_t charge;
  Double_t trkPt;

  TkMuon()
  {
    Init();
  }

  TkMuon(NtupleHandle* ntuple, Int_t index): TkMuon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt   = ntuple->TkMuon_pt[index];
    eta  = ntuple->TkMuon_eta[index];
    phi  = ntuple->TkMuon_phi[index];
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    charge = ntuple->TkMuon_charge[index];
    trkPt  = ntuple->TkMuon_trkPt[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    charge = -999;
    trkPt = -999;
  }
};

class L1Muon : public Object
{
public:
  Double_t charge;
  Double_t quality;

  L1Muon()
  {
    Init();
  }

  L1Muon(NtupleHandle* ntuple, Int_t index): L1Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt   = ntuple->L1Muon_pt[index];
    eta  = ntuple->L1Muon_eta[index];
    phi  = ntuple->L1Muon_phi[index];
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    charge  = ntuple->L1Muon_charge[index];
    quality = ntuple->L1Muon_quality[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    charge = -999;
    quality = -999;
  }
};

class IterL3OIMuon : public Object
{
public:
  Double_t inner_pt;
  Double_t inner_eta;
  Double_t inner_phi;
  Double_t inner_charge;

  Double_t outer_pt;
  Double_t outer_eta;
  Double_t outer_phi;
  Double_t outer_charge;

  Double_t global_pt;
  Double_t global_eta;
  Double_t global_phi;
  Double_t global_charge;

  IterL3OIMuon()
  {
    Init();
  }

  IterL3OIMuon(NtupleHandle* ntuple, Int_t index): IterL3OIMuon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    inner_pt     = ntuple->iterL3OI_inner_pt[index];
    inner_eta    = ntuple->iterL3OI_inner_eta[index];
    inner_phi    = ntuple->iterL3OI_inner_phi[index];
    inner_charge = ntuple->iterL3OI_inner_charge[index];

    pt = inner_pt;
    eta = inner_eta;
    phi = inner_phi;
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    outer_pt     = ntuple->iterL3OI_outer_pt[index];
    outer_eta    = ntuple->iterL3OI_outer_eta[index];
    outer_phi    = ntuple->iterL3OI_outer_phi[index];
    outer_charge = ntuple->iterL3OI_outer_charge[index];

    global_pt     = ntuple->iterL3OI_global_pt[index];
    global_eta    = ntuple->iterL3OI_global_eta[index];
    global_phi    = ntuple->iterL3OI_global_phi[index];
    global_charge = ntuple->iterL3OI_global_charge[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;

    inner_pt = -999;
    inner_eta = -999;
    inner_phi = -999;
    inner_charge = -999;

    outer_pt = -999;
    outer_eta = -999;
    outer_phi = -999;
    outer_charge = -999;

    global_pt = -999;
    global_eta = -999;
    global_phi = -999;
    global_charge = -999;
  }
};

class IterL3IOFromL2Muon : public Object
{
public:
  Double_t inner_pt;
  Double_t inner_eta;
  Double_t inner_phi;
  Double_t inner_charge;

  Double_t outer_pt;
  Double_t outer_eta;
  Double_t outer_phi;
  Double_t outer_charge;

  Double_t global_pt;
  Double_t global_eta;
  Double_t global_phi;
  Double_t global_charge;

  IterL3IOFromL2Muon()
  {
    Init();
  }

  IterL3IOFromL2Muon(NtupleHandle* ntuple, Int_t index): IterL3IOFromL2Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    inner_pt     = ntuple->iterL3IOFromL2_inner_pt[index];
    inner_eta    = ntuple->iterL3IOFromL2_inner_eta[index];
    inner_phi    = ntuple->iterL3IOFromL2_inner_phi[index];
    inner_charge = ntuple->iterL3IOFromL2_inner_charge[index];

    pt = inner_pt;
    eta = inner_eta;
    phi = inner_phi;
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    outer_pt     = ntuple->iterL3IOFromL2_outer_pt[index];
    outer_eta    = ntuple->iterL3IOFromL2_outer_eta[index];
    outer_phi    = ntuple->iterL3IOFromL2_outer_phi[index];
    outer_charge = ntuple->iterL3IOFromL2_outer_charge[index];

    global_pt     = ntuple->iterL3IOFromL2_global_pt[index];
    global_eta    = ntuple->iterL3IOFromL2_global_eta[index];
    global_phi    = ntuple->iterL3IOFromL2_global_phi[index];
    global_charge = ntuple->iterL3IOFromL2_global_charge[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;

    inner_pt = -999;
    inner_eta = -999;
    inner_phi = -999;
    inner_charge = -999;

    outer_pt = -999;
    outer_eta = -999;
    outer_phi = -999;
    outer_charge = -999;

    global_pt = -999;
    global_eta = -999;
    global_phi = -999;
    global_charge = -999;
  }
};

class IterL3FromL2Muon : public Object
{
public:
  Double_t inner_pt;
  Double_t inner_eta;
  Double_t inner_phi;
  Double_t inner_charge;

  Double_t outer_pt;
  Double_t outer_eta;
  Double_t outer_phi;
  Double_t outer_charge;

  Double_t global_pt;
  Double_t global_eta;
  Double_t global_phi;
  Double_t global_charge;

  IterL3FromL2Muon()
  {
    Init();
  }

  IterL3FromL2Muon(NtupleHandle* ntuple, Int_t index): IterL3FromL2Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    inner_pt     = ntuple->iterL3FromL2_inner_pt[index];
    inner_eta    = ntuple->iterL3FromL2_inner_eta[index];
    inner_phi    = ntuple->iterL3FromL2_inner_phi[index];
    inner_charge = ntuple->iterL3FromL2_inner_charge[index];

    pt = inner_pt;
    eta = inner_eta;
    phi = inner_phi;
    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    outer_pt     = ntuple->iterL3FromL2_outer_pt[index];
    outer_eta    = ntuple->iterL3FromL2_outer_eta[index];
    outer_phi    = ntuple->iterL3FromL2_outer_phi[index];
    outer_charge = ntuple->iterL3FromL2_outer_charge[index];

    global_pt     = ntuple->iterL3FromL2_global_pt[index];
    global_eta    = ntuple->iterL3FromL2_global_eta[index];
    global_phi    = ntuple->iterL3FromL2_global_phi[index];
    global_charge = ntuple->iterL3FromL2_global_charge[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;

    inner_pt = -999;
    inner_eta = -999;
    inner_phi = -999;
    inner_charge = -999;

    outer_pt = -999;
    outer_eta = -999;
    outer_phi = -999;
    outer_charge = -999;

    global_pt = -999;
    global_eta = -999;
    global_phi = -999;
    global_charge = -999;
  }
};

class IterL3IOFromL1Muon : public Object
{
public:
  Double_t charge;

  IterL3IOFromL1Muon()
  {
    Init();
  }

  IterL3IOFromL1Muon(NtupleHandle* ntuple, Int_t index): IterL3IOFromL1Muon()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt     = ntuple->iterL3IOFromL1_pt[index];
    eta    = ntuple->iterL3IOFromL1_eta[index];
    phi    = ntuple->iterL3IOFromL1_phi[index];
    charge = ntuple->iterL3IOFromL1_charge[index];

    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);
  }

private:
  void Init()
  {
    pt = 0;
    eta = 0;
    phi = 0;
    charge = 0;
    mass = 0;
  }
};

class IterL3MuonNoID : public Object
{
public:
  Double_t pt;
  Double_t eta;
  Double_t phi;
  Double_t charge;
  Int_t    isGLB;
  Int_t    isSTA;
  Int_t    isTRK;

  IterL3MuonNoID()
  {
    Init();
  }

  IterL3MuonNoID(NtupleHandle* ntuple, Int_t index): IterL3MuonNoID()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    pt     = ntuple->iterL3MuonNoID_pt[index];
    eta    = ntuple->iterL3MuonNoID_eta[index];
    phi    = ntuple->iterL3MuonNoID_phi[index];
    charge = ntuple->iterL3MuonNoID_charge[index];

    mass = MuonHLT::M_mu;
    vecP.SetPtEtaPhiM(pt, eta, phi, mass);

    isGLB = ntuple->iterL3MuonNoID_isGLB[index];
    isSTA = ntuple->iterL3MuonNoID_isSTA[index];
    isTRK = ntuple->iterL3MuonNoID_isTRK[index];
  }

private:
  void Init()
  {
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
    charge = -999;
    isGLB = -999;
    isSTA = -999;
    isTRK = -999;
  }
};


class HLTObject : public Object
{
public:
  TString filterName;

  HLTObject()
  {
    Init();
  }

  HLTObject(NtupleHandle* ntuple, Int_t index): HLTObject()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    filterName = ntuple->vec_filterName->at((UInt_t)index);
    pt  = ntuple->vec_HLTObj_pt->at((UInt_t)index);
    eta = ntuple->vec_HLTObj_eta->at((UInt_t)index);
    phi = ntuple->vec_HLTObj_phi->at((UInt_t)index);
    mass = MuonHLT::M_mu; // -- it could be wrong if you don't use muon HLT

    vecP.SetPtEtaPhiM(pt, eta, phi, mass);
  }

private:
  void Init()
  {
    filterName = "";
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
  }
};

class MYHLTObject : public Object
{
public:
  TString filterName;

  MYHLTObject()
  {
    Init();
  }

  MYHLTObject(NtupleHandle* ntuple, Int_t index): MYHLTObject()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {
    filterName = ntuple->vec_myFilterName->at((UInt_t)index);
    pt  = ntuple->vec_myHLTObj_pt->at((UInt_t)index);
    eta = ntuple->vec_myHLTObj_eta->at((UInt_t)index);
    phi = ntuple->vec_myHLTObj_phi->at((UInt_t)index);
    mass = MuonHLT::M_mu; // -- it could be wrong if you don't use muon HLT

    vecP.SetPtEtaPhiM(pt, eta, phi, mass);
  }

private:
  void Init()
  {
    filterName = "";
    pt = -999;
    eta = -999;
    phi = -999;
    mass = -999;
  }
};

// -- template
class ParticleTemplate : public Object
{
public:
  ParticleTemplate()
  {
    Init();
  }

  ParticleTemplate(NtupleHandle* ntuple, Int_t index): ParticleTemplate()
  {
    FillVariable(ntuple, index);
  }

  void FillVariable(NtupleHandle* ntuple, Int_t index)
  {

  }

private:
  void Init()
  {

  }
};

}; // -- end of MuonHLT namespace