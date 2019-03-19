#pragma once

#define ArrSize 50000
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <vector>

namespace MuonHLT
{

// -- for the physical variables, do not use _ at the end of the variable name.
// -- because getter functions for those variables will not be created (too many!)
class NtupleHandle
{
public:
  TChain *chain_;

  // -- event information
  Bool_t          isRealData;
  Int_t           runNum;
  Int_t           lumiBlockNum;
  ULong64_t       eventNum;
  Int_t           nVertex;
  Double_t        bunchID;
  Double_t        instLumi;
  Double_t        dataPU;
  Double_t        dataPURMS;
  Double_t        bunchLumi;
  Double_t        offlineInstLumi;
  Double_t        offlineDataPU;
  Double_t        offlineDataPURMS;
  Double_t        offlineBunchLumi;
  Int_t           truePU;
  Double_t        genEventWeight;

  // -- generator inforomation
  Int_t           nGenParticle;
  Int_t           genParticle_ID[ArrSize];
  Int_t           genParticle_status[ArrSize];
  Int_t           genParticle_mother[ArrSize];
  Double_t        genParticle_pt[ArrSize];
  Double_t        genParticle_eta[ArrSize];
  Double_t        genParticle_phi[ArrSize];
  Double_t        genParticle_px[ArrSize];
  Double_t        genParticle_py[ArrSize];
  Double_t        genParticle_pz[ArrSize];
  Double_t        genParticle_energy[ArrSize];
  Double_t        genParticle_charge[ArrSize];
  Int_t           genParticle_isPrompt[ArrSize];
  Int_t           genParticle_isPromptFinalState[ArrSize];
  Int_t           genParticle_isTauDecayProduct[ArrSize];
  Int_t           genParticle_isPromptTauDecayProduct[ArrSize];
  Int_t           genParticle_isDirectPromptTauDecayProductFinalState[ArrSize];
  Int_t           genParticle_isHardProcess[ArrSize];
  Int_t           genParticle_isLastCopy[ArrSize];
  Int_t           genParticle_isLastCopyBeforeFSR[ArrSize];
  Int_t           genParticle_isPromptDecayed[ArrSize];
  Int_t           genParticle_isDecayedLeptonHadron[ArrSize];
  Int_t           genParticle_fromHardProcessBeforeFSR[ArrSize];
  Int_t           genParticle_fromHardProcessDecayed[ArrSize];
  Int_t           genParticle_fromHardProcessFinalState[ArrSize];
  Int_t           genParticle_isMostlyLikePythia6Status3[ArrSize];

  // -- trigger information
  vector<string>  *vec_firedTrigger;
  vector<string>  *vec_filterName;
  vector<double>  *vec_HLTObj_pt;
  vector<double>  *vec_HLTObj_eta;
  vector<double>  *vec_HLTObj_phi;
  vector<string>  *vec_myFiredTrigger;
  vector<string>  *vec_myFilterName;
  vector<double>  *vec_myHLTObj_pt;
  vector<double>  *vec_myHLTObj_eta;
  vector<double>  *vec_myHLTObj_phi;

  // -- offline muons
  Int_t           nMuon;
  Double_t        muon_pt[ArrSize];
  Double_t        muon_eta[ArrSize];
  Double_t        muon_phi[ArrSize];
  Double_t        muon_px[ArrSize];
  Double_t        muon_py[ArrSize];
  Double_t        muon_pz[ArrSize];
  Double_t        muon_dB[ArrSize];
  Double_t        muon_charge[ArrSize];
  Int_t           muon_isGLB[ArrSize];
  Int_t           muon_isSTA[ArrSize];
  Int_t           muon_isTRK[ArrSize];
  Int_t           muon_isPF[ArrSize];
  Int_t           muon_isTight[ArrSize];
  Int_t           muon_isMedium[ArrSize];
  Int_t           muon_isLoose[ArrSize];
  Int_t           muon_isHighPt[ArrSize];
  Int_t           muon_isSoft[ArrSize];
  Double_t        muon_iso03_sumPt[ArrSize];
  Double_t        muon_iso03_hadEt[ArrSize];
  Double_t        muon_iso03_emEt[ArrSize];
  Double_t        muon_PFIso03_charged[ArrSize];
  Double_t        muon_PFIso03_neutral[ArrSize];
  Double_t        muon_PFIso03_photon[ArrSize];
  Double_t        muon_PFIso03_sumPU[ArrSize];
  Double_t        muon_PFIso04_charged[ArrSize];
  Double_t        muon_PFIso04_neutral[ArrSize];
  Double_t        muon_PFIso04_photon[ArrSize];
  Double_t        muon_PFIso04_sumPU[ArrSize];
  Double_t        muon_PFCluster03_ECAL[ArrSize];
  Double_t        muon_PFCluster03_HCAL[ArrSize];
  Double_t        muon_PFCluster04_ECAL[ArrSize];
  Double_t        muon_PFCluster04_HCAL[ArrSize];
  Double_t        muon_normChi2_global[ArrSize];
  Int_t           muon_nTrackerHit_global[ArrSize];
  Int_t           muon_nTrackerLayer_global[ArrSize];
  Int_t           muon_nPixelHit_global[ArrSize];
  Int_t           muon_nMuonHit_global[ArrSize];
  Double_t        muon_normChi2_inner[ArrSize];
  Int_t           muon_nTrackerHit_inner[ArrSize];
  Int_t           muon_nTrackerLayer_inner[ArrSize];
  Int_t           muon_nPixelHit_inner[ArrSize];
  Double_t        muon_pt_tuneP[ArrSize];
  Double_t        muon_ptError_tuneP[ArrSize];
  Double_t        muon_dxyVTX_best[ArrSize];
  Double_t        muon_dzVTX_best[ArrSize];
  Int_t           muon_nMatchedStation[ArrSize];
  Int_t           muon_nMatchedRPCLayer[ArrSize];
  Int_t           muon_stationMask[ArrSize];

  // -- L3 muons
  Int_t           nL3Muon;
  Double_t        L3Muon_pt[ArrSize];
  Double_t        L3Muon_eta[ArrSize];
  Double_t        L3Muon_phi[ArrSize];
  Double_t        L3Muon_charge[ArrSize];
  Double_t        L3Muon_trkPt[ArrSize];

  // -- L2 muons
  Int_t           nL2Muon;
  Double_t        L2Muon_pt[ArrSize];
  Double_t        L2Muon_eta[ArrSize];
  Double_t        L2Muon_phi[ArrSize];
  Double_t        L2Muon_charge[ArrSize];
  Double_t        L2Muon_trkPt[ArrSize];

  // -- Tk muons
  Int_t           nTkMuon;
  Double_t        TkMuon_pt[ArrSize];
  Double_t        TkMuon_eta[ArrSize];
  Double_t        TkMuon_phi[ArrSize];
  Double_t        TkMuon_charge[ArrSize];
  Double_t        TkMuon_trkPt[ArrSize];

  // -- L1 muons
  Int_t           nL1Muon;
  Double_t        L1Muon_pt[ArrSize];
  Double_t        L1Muon_eta[ArrSize];
  Double_t        L1Muon_phi[ArrSize];
  Double_t        L1Muon_charge[ArrSize];
  Double_t        L1Muon_quality[ArrSize];

  // -- iterative L3: outside-in (from L2)
  Int_t           nIterL3OI;
  Double_t        iterL3OI_inner_pt[ArrSize];
  Double_t        iterL3OI_inner_eta[ArrSize];
  Double_t        iterL3OI_inner_phi[ArrSize];
  Double_t        iterL3OI_inner_charge[ArrSize];
  Double_t        iterL3OI_outer_pt[ArrSize];
  Double_t        iterL3OI_outer_eta[ArrSize];
  Double_t        iterL3OI_outer_phi[ArrSize];
  Double_t        iterL3OI_outer_charge[ArrSize];
  Double_t        iterL3OI_global_pt[ArrSize];
  Double_t        iterL3OI_global_eta[ArrSize];
  Double_t        iterL3OI_global_phi[ArrSize];
  Double_t        iterL3OI_global_charge[ArrSize];

  // -- iterative L3: inside-out from L2
  Int_t           nIterL3IOFromL2;
  Double_t        iterL3IOFromL2_inner_pt[ArrSize];
  Double_t        iterL3IOFromL2_inner_eta[ArrSize];
  Double_t        iterL3IOFromL2_inner_phi[ArrSize];
  Double_t        iterL3IOFromL2_inner_charge[ArrSize];
  Double_t        iterL3IOFromL2_outer_pt[ArrSize];
  Double_t        iterL3IOFromL2_outer_eta[ArrSize];
  Double_t        iterL3IOFromL2_outer_phi[ArrSize];
  Double_t        iterL3IOFromL2_outer_charge[ArrSize];
  Double_t        iterL3IOFromL2_global_pt[ArrSize];
  Double_t        iterL3IOFromL2_global_eta[ArrSize];
  Double_t        iterL3IOFromL2_global_phi[ArrSize];
  Double_t        iterL3IOFromL2_global_charge[ArrSize];

  // -- iterative L3: outside-in + inside-out from L2
  Int_t           nIterL3FromL2;
  Double_t        iterL3FromL2_inner_pt[ArrSize];
  Double_t        iterL3FromL2_inner_eta[ArrSize];
  Double_t        iterL3FromL2_inner_phi[ArrSize];
  Double_t        iterL3FromL2_inner_charge[ArrSize];
  Double_t        iterL3FromL2_outer_pt[ArrSize];
  Double_t        iterL3FromL2_outer_eta[ArrSize];
  Double_t        iterL3FromL2_outer_phi[ArrSize];
  Double_t        iterL3FromL2_outer_charge[ArrSize];
  Double_t        iterL3FromL2_global_pt[ArrSize];
  Double_t        iterL3FromL2_global_eta[ArrSize];
  Double_t        iterL3FromL2_global_phi[ArrSize];
  Double_t        iterL3FromL2_global_charge[ArrSize];

  // -- iterative L3: inside-out from L1
  Int_t           nIterL3IOFromL1;
  Double_t        iterL3IOFromL1_pt[ArrSize];
  Double_t        iterL3IOFromL1_eta[ArrSize];
  Double_t        iterL3IOFromL1_phi[ArrSize];
  Double_t        iterL3IOFromL1_charge[ArrSize];

  // -- iterative L3: muon before ID filter
  Int_t           nIterL3MuonNoID;
  Double_t        iterL3MuonNoID_pt[ArrSize];
  Double_t        iterL3MuonNoID_eta[ArrSize];
  Double_t        iterL3MuonNoID_phi[ArrSize];
  Double_t        iterL3MuonNoID_charge[ArrSize];
  Int_t           iterL3MuonNoID_isGLB[ArrSize];
  Int_t           iterL3MuonNoID_isSTA[ArrSize];
  Int_t           iterL3MuonNoID_isTRK[ArrSize];


  NtupleHandle()
  {
    // -- init. vectors
    vec_firedTrigger = 0;
    vec_filterName = 0;
    vec_HLTObj_pt = 0;
    vec_HLTObj_eta = 0;
    vec_HLTObj_phi = 0;
    vec_myFiredTrigger = 0;
    vec_myFilterName = 0;
    vec_myHLTObj_pt = 0;
    vec_myHLTObj_eta = 0;
    vec_myHLTObj_phi = 0;
  }

  NtupleHandle(TChain* chain): NtupleHandle()
  {
    chain_ = chain;
    chain_->SetBranchStatus("*", 0);

    TurnOnBranches_Event();
    TurnOnBranches_Trigger();
    // TurnOnBranches_GenParticle();
    // TurnOnBranches_Muon();
    // TurnOnBranches_HLTMuon();
    // TurnOnBranches_IterL3Muon();
  }

  void GetEvent(Int_t index)
  {
    chain_->GetEntry(index);
  }

  void TurnOnBranches_Event()
  {
    chain_->SetBranchStatus("isRealData", 1);
    chain_->SetBranchStatus("runNum", 1);
    chain_->SetBranchStatus("lumiBlockNum", 1);
    chain_->SetBranchStatus("eventNum", 1);
    chain_->SetBranchStatus("nVertex", 1);
    chain_->SetBranchStatus("bunchID", 1);
    chain_->SetBranchStatus("instLumi", 1);
    chain_->SetBranchStatus("dataPU", 1);
    chain_->SetBranchStatus("dataPURMS", 1);
    chain_->SetBranchStatus("bunchLumi", 1);
    chain_->SetBranchStatus("offlineInstLumi", 1);
    chain_->SetBranchStatus("offlineDataPU", 1);
    chain_->SetBranchStatus("offlineDataPURMS", 1);
    chain_->SetBranchStatus("offlineBunchLumi", 1);
    chain_->SetBranchStatus("truePU", 1);
    chain_->SetBranchStatus("genEventWeight", 1);

    chain_->SetBranchAddress("isRealData", &isRealData);
    chain_->SetBranchAddress("runNum", &runNum);
    chain_->SetBranchAddress("lumiBlockNum", &lumiBlockNum);
    chain_->SetBranchAddress("eventNum", &eventNum);
    chain_->SetBranchAddress("nVertex", &nVertex);
    chain_->SetBranchAddress("bunchID", &bunchID);
    chain_->SetBranchAddress("instLumi", &instLumi);
    chain_->SetBranchAddress("dataPU", &dataPU);
    chain_->SetBranchAddress("dataPURMS", &dataPURMS);
    chain_->SetBranchAddress("bunchLumi", &bunchLumi);
    chain_->SetBranchAddress("offlineInstLumi", &offlineInstLumi);
    chain_->SetBranchAddress("offlineDataPU", &offlineDataPU);
    chain_->SetBranchAddress("offlineDataPURMS", &offlineDataPURMS);
    chain_->SetBranchAddress("offlineBunchLumi", &offlineBunchLumi);
    chain_->SetBranchAddress("truePU", &truePU);
    chain_->SetBranchAddress("genEventWeight", &genEventWeight);
  }

  void TurnOnBranches_Trigger()
  {
    chain_->SetBranchStatus("vec_firedTrigger", 1);
    chain_->SetBranchAddress("vec_firedTrigger", &vec_firedTrigger);

    chain_->SetBranchStatus("vec_filterName", 1);
    chain_->SetBranchAddress("vec_filterName", &vec_filterName);

    chain_->SetBranchStatus("vec_HLTObj_pt", 1);
    chain_->SetBranchAddress("vec_HLTObj_pt", &vec_HLTObj_pt);

    chain_->SetBranchStatus("vec_HLTObj_eta", 1);
    chain_->SetBranchAddress("vec_HLTObj_eta", &vec_HLTObj_eta);

    chain_->SetBranchStatus("vec_HLTObj_phi", 1);
    chain_->SetBranchAddress("vec_HLTObj_phi", &vec_HLTObj_phi);

    chain_->SetBranchStatus("vec_myFiredTrigger", 1);
    chain_->SetBranchAddress("vec_myFiredTrigger", &vec_myFiredTrigger);

    chain_->SetBranchStatus("vec_myFilterName", 1);
    chain_->SetBranchAddress("vec_myFilterName", &vec_myFilterName);

    chain_->SetBranchStatus("vec_myHLTObj_pt", 1);
    chain_->SetBranchAddress("vec_myHLTObj_pt", &vec_myHLTObj_pt);

    chain_->SetBranchStatus("vec_myHLTObj_eta", 1);
    chain_->SetBranchAddress("vec_myHLTObj_eta", &vec_myHLTObj_eta);

    chain_->SetBranchStatus("vec_myHLTObj_phi", 1);
    chain_->SetBranchAddress("vec_myHLTObj_phi", &vec_myHLTObj_phi);
  }

  void TurnOnBranches_GenParticle()
  {
    chain_->SetBranchStatus("nGenParticle", 1);
    chain_->SetBranchAddress("nGenParticle", &nGenParticle);

    chain_->SetBranchStatus("genParticle_ID", 1);
    chain_->SetBranchAddress("genParticle_ID", &genParticle_ID);

    chain_->SetBranchStatus("genParticle_status", 1);
    chain_->SetBranchAddress("genParticle_status", &genParticle_status);

    chain_->SetBranchStatus("genParticle_mother", 1);
    chain_->SetBranchAddress("genParticle_mother", &genParticle_mother);

    chain_->SetBranchStatus("genParticle_pt", 1);
    chain_->SetBranchAddress("genParticle_pt", &genParticle_pt);

    chain_->SetBranchStatus("genParticle_eta", 1);
    chain_->SetBranchAddress("genParticle_eta", &genParticle_eta);

    chain_->SetBranchStatus("genParticle_phi", 1);
    chain_->SetBranchAddress("genParticle_phi", &genParticle_phi);

    chain_->SetBranchStatus("genParticle_px", 1);
    chain_->SetBranchAddress("genParticle_px", &genParticle_px);

    chain_->SetBranchStatus("genParticle_py", 1);
    chain_->SetBranchAddress("genParticle_py", &genParticle_py);

    chain_->SetBranchStatus("genParticle_pz", 1);
    chain_->SetBranchAddress("genParticle_pz", &genParticle_pz);

    chain_->SetBranchStatus("genParticle_energy", 1);
    chain_->SetBranchAddress("genParticle_energy", &genParticle_energy);

    chain_->SetBranchStatus("genParticle_charge", 1);
    chain_->SetBranchAddress("genParticle_charge", &genParticle_charge);

    chain_->SetBranchStatus("genParticle_isPrompt", 1);
    chain_->SetBranchAddress("genParticle_isPrompt", &genParticle_isPrompt);

    chain_->SetBranchStatus("genParticle_isPromptFinalState", 1);
    chain_->SetBranchAddress("genParticle_isPromptFinalState", &genParticle_isPromptFinalState);

    chain_->SetBranchStatus("genParticle_isTauDecayProduct", 1);
    chain_->SetBranchAddress("genParticle_isTauDecayProduct", &genParticle_isTauDecayProduct);

    chain_->SetBranchStatus("genParticle_isPromptTauDecayProduct", 1);
    chain_->SetBranchAddress("genParticle_isPromptTauDecayProduct", &genParticle_isPromptTauDecayProduct);

    chain_->SetBranchStatus("genParticle_isDirectPromptTauDecayProductFinalState", 1);
    chain_->SetBranchAddress("genParticle_isDirectPromptTauDecayProductFinalState", &genParticle_isDirectPromptTauDecayProductFinalState);

    chain_->SetBranchStatus("genParticle_isHardProcess", 1);
    chain_->SetBranchAddress("genParticle_isHardProcess", &genParticle_isHardProcess);

    chain_->SetBranchStatus("genParticle_isLastCopy", 1);
    chain_->SetBranchAddress("genParticle_isLastCopy", &genParticle_isLastCopy);

    chain_->SetBranchStatus("genParticle_isLastCopyBeforeFSR", 1);
    chain_->SetBranchAddress("genParticle_isLastCopyBeforeFSR", &genParticle_isLastCopyBeforeFSR);

    chain_->SetBranchStatus("genParticle_isPromptDecayed", 1);
    chain_->SetBranchAddress("genParticle_isPromptDecayed", &genParticle_isPromptDecayed);

    chain_->SetBranchStatus("genParticle_isDecayedLeptonHadron", 1);
    chain_->SetBranchAddress("genParticle_isDecayedLeptonHadron", &genParticle_isDecayedLeptonHadron);

    chain_->SetBranchStatus("genParticle_fromHardProcessBeforeFSR", 1);
    chain_->SetBranchAddress("genParticle_fromHardProcessBeforeFSR", &genParticle_fromHardProcessBeforeFSR);

    chain_->SetBranchStatus("genParticle_fromHardProcessDecayed", 1);
    chain_->SetBranchAddress("genParticle_fromHardProcessDecayed", &genParticle_fromHardProcessDecayed);

    chain_->SetBranchStatus("genParticle_fromHardProcessFinalState", 1);
    chain_->SetBranchAddress("genParticle_fromHardProcessFinalState", &genParticle_fromHardProcessFinalState);

    chain_->SetBranchStatus("genParticle_isMostlyLikePythia6Status3", 1);
    chain_->SetBranchAddress("genParticle_isMostlyLikePythia6Status3", &genParticle_isMostlyLikePythia6Status3);
  }

  void TurnOnBranches_Muon()
  {
    chain_->SetBranchStatus("nMuon", 1);
    chain_->SetBranchAddress("nMuon", &nMuon);

    chain_->SetBranchStatus("muon_pt", 1);
    chain_->SetBranchAddress("muon_pt", &muon_pt);

    chain_->SetBranchStatus("muon_eta", 1);
    chain_->SetBranchAddress("muon_eta", &muon_eta);

    chain_->SetBranchStatus("muon_phi", 1);
    chain_->SetBranchAddress("muon_phi", &muon_phi);

    chain_->SetBranchStatus("muon_px", 1);
    chain_->SetBranchAddress("muon_px", &muon_px);

    chain_->SetBranchStatus("muon_py", 1);
    chain_->SetBranchAddress("muon_py", &muon_py);

    chain_->SetBranchStatus("muon_pz", 1);
    chain_->SetBranchAddress("muon_pz", &muon_pz);

    chain_->SetBranchStatus("muon_dB", 1);
    chain_->SetBranchAddress("muon_dB", &muon_dB);

    chain_->SetBranchStatus("muon_charge", 1);
    chain_->SetBranchAddress("muon_charge", &muon_charge);

    chain_->SetBranchStatus("muon_isGLB", 1);
    chain_->SetBranchAddress("muon_isGLB", &muon_isGLB);

    chain_->SetBranchStatus("muon_isSTA", 1);
    chain_->SetBranchAddress("muon_isSTA", &muon_isSTA);

    chain_->SetBranchStatus("muon_isTRK", 1);
    chain_->SetBranchAddress("muon_isTRK", &muon_isTRK);

    chain_->SetBranchStatus("muon_isPF", 1);
    chain_->SetBranchAddress("muon_isPF", &muon_isPF);

    chain_->SetBranchStatus("muon_isTight", 1);
    chain_->SetBranchAddress("muon_isTight", &muon_isTight);

    chain_->SetBranchStatus("muon_isMedium", 1);
    chain_->SetBranchAddress("muon_isMedium", &muon_isMedium);

    chain_->SetBranchStatus("muon_isLoose", 1);
    chain_->SetBranchAddress("muon_isLoose", &muon_isLoose);

    chain_->SetBranchStatus("muon_isHighPt", 1);
    chain_->SetBranchAddress("muon_isHighPt", &muon_isHighPt);

    chain_->SetBranchStatus("muon_isSoft", 1);
    chain_->SetBranchAddress("muon_isSoft", &muon_isSoft);

    chain_->SetBranchStatus("muon_iso03_sumPt", 1);
    chain_->SetBranchAddress("muon_iso03_sumPt", &muon_iso03_sumPt);

    chain_->SetBranchStatus("muon_iso03_hadEt", 1);
    chain_->SetBranchAddress("muon_iso03_hadEt", &muon_iso03_hadEt);

    chain_->SetBranchStatus("muon_iso03_emEt", 1);
    chain_->SetBranchAddress("muon_iso03_emEt", &muon_iso03_emEt);

    chain_->SetBranchStatus("muon_PFIso03_charged", 1);
    chain_->SetBranchAddress("muon_PFIso03_charged", &muon_PFIso03_charged);

    chain_->SetBranchStatus("muon_PFIso03_neutral", 1);
    chain_->SetBranchAddress("muon_PFIso03_neutral", &muon_PFIso03_neutral);

    chain_->SetBranchStatus("muon_PFIso03_photon", 1);
    chain_->SetBranchAddress("muon_PFIso03_photon", &muon_PFIso03_photon);

    chain_->SetBranchStatus("muon_PFIso03_sumPU", 1);
    chain_->SetBranchAddress("muon_PFIso03_sumPU", &muon_PFIso03_sumPU);

    chain_->SetBranchStatus("muon_PFIso04_charged", 1);
    chain_->SetBranchAddress("muon_PFIso04_charged", &muon_PFIso04_charged);

    chain_->SetBranchStatus("muon_PFIso04_neutral", 1);
    chain_->SetBranchAddress("muon_PFIso04_neutral", &muon_PFIso04_neutral);

    chain_->SetBranchStatus("muon_PFIso04_photon", 1);
    chain_->SetBranchAddress("muon_PFIso04_photon", &muon_PFIso04_photon);

    chain_->SetBranchStatus("muon_PFIso04_sumPU", 1);
    chain_->SetBranchAddress("muon_PFIso04_sumPU", &muon_PFIso04_sumPU);

    chain_->SetBranchStatus("muon_PFCluster03_ECAL", 1);
    chain_->SetBranchAddress("muon_PFCluster03_ECAL", &muon_PFCluster03_ECAL);

    chain_->SetBranchStatus("muon_PFCluster03_HCAL", 1);
    chain_->SetBranchAddress("muon_PFCluster03_HCAL", &muon_PFCluster03_HCAL);

    chain_->SetBranchStatus("muon_PFCluster04_ECAL", 1);
    chain_->SetBranchAddress("muon_PFCluster04_ECAL", &muon_PFCluster04_ECAL);

    chain_->SetBranchStatus("muon_PFCluster04_HCAL", 1);
    chain_->SetBranchAddress("muon_PFCluster04_HCAL", &muon_PFCluster04_HCAL);

    chain_->SetBranchStatus("muon_normChi2_global", 1);
    chain_->SetBranchAddress("muon_normChi2_global", &muon_normChi2_global);

    chain_->SetBranchStatus("muon_nTrackerHit_global", 1);
    chain_->SetBranchAddress("muon_nTrackerHit_global", &muon_nTrackerHit_global);

    chain_->SetBranchStatus("muon_nTrackerLayer_global", 1);
    chain_->SetBranchAddress("muon_nTrackerLayer_global", &muon_nTrackerLayer_global);

    chain_->SetBranchStatus("muon_nPixelHit_global", 1);
    chain_->SetBranchAddress("muon_nPixelHit_global", &muon_nPixelHit_global);

    chain_->SetBranchStatus("muon_nMuonHit_global", 1);
    chain_->SetBranchAddress("muon_nMuonHit_global", &muon_nMuonHit_global);

    chain_->SetBranchStatus("muon_normChi2_inner", 1);
    chain_->SetBranchAddress("muon_normChi2_inner", &muon_normChi2_inner);

    chain_->SetBranchStatus("muon_nTrackerHit_inner", 1);
    chain_->SetBranchAddress("muon_nTrackerHit_inner", &muon_nTrackerHit_inner);

    chain_->SetBranchStatus("muon_nTrackerLayer_inner", 1);
    chain_->SetBranchAddress("muon_nTrackerLayer_inner", &muon_nTrackerLayer_inner);

    chain_->SetBranchStatus("muon_nPixelHit_inner", 1);
    chain_->SetBranchAddress("muon_nPixelHit_inner", &muon_nPixelHit_inner);

    chain_->SetBranchStatus("muon_pt_tuneP", 1);
    chain_->SetBranchAddress("muon_pt_tuneP", &muon_pt_tuneP);

    chain_->SetBranchStatus("muon_ptError_tuneP", 1);
    chain_->SetBranchAddress("muon_ptError_tuneP", &muon_ptError_tuneP);

    chain_->SetBranchStatus("muon_dxyVTX_best", 1);
    chain_->SetBranchAddress("muon_dxyVTX_best", &muon_dxyVTX_best);

    chain_->SetBranchStatus("muon_dzVTX_best", 1);
    chain_->SetBranchAddress("muon_dzVTX_best", &muon_dzVTX_best);

    chain_->SetBranchStatus("muon_nMatchedStation", 1);
    chain_->SetBranchAddress("muon_nMatchedStation", &muon_nMatchedStation);

    chain_->SetBranchStatus("muon_nMatchedRPCLayer", 1);
    chain_->SetBranchAddress("muon_nMatchedRPCLayer", &muon_nMatchedRPCLayer);

    chain_->SetBranchStatus("muon_stationMask", 1);
    chain_->SetBranchAddress("muon_stationMask", &muon_stationMask);
  }

  void TurnOnBranches_HLTMuon()
  {
    chain_->SetBranchStatus("nL3Muon", 1);
    chain_->SetBranchAddress("nL3Muon", &nL3Muon);

    chain_->SetBranchStatus("L3Muon_pt", 1);
    chain_->SetBranchAddress("L3Muon_pt", &L3Muon_pt);

    chain_->SetBranchStatus("L3Muon_eta", 1);
    chain_->SetBranchAddress("L3Muon_eta", &L3Muon_eta);

    chain_->SetBranchStatus("L3Muon_phi", 1);
    chain_->SetBranchAddress("L3Muon_phi", &L3Muon_phi);

    chain_->SetBranchStatus("L3Muon_charge", 1);
    chain_->SetBranchAddress("L3Muon_charge", &L3Muon_charge);

    chain_->SetBranchStatus("L3Muon_trkPt", 1);
    chain_->SetBranchAddress("L3Muon_trkPt", &L3Muon_trkPt);

    chain_->SetBranchStatus("nL2Muon", 1);
    chain_->SetBranchAddress("nL2Muon", &nL2Muon);

    chain_->SetBranchStatus("L2Muon_pt", 1);
    chain_->SetBranchAddress("L2Muon_pt", &L2Muon_pt);

    chain_->SetBranchStatus("L2Muon_eta", 1);
    chain_->SetBranchAddress("L2Muon_eta", &L2Muon_eta);

    chain_->SetBranchStatus("L2Muon_phi", 1);
    chain_->SetBranchAddress("L2Muon_phi", &L2Muon_phi);

    chain_->SetBranchStatus("L2Muon_charge", 1);
    chain_->SetBranchAddress("L2Muon_charge", &L2Muon_charge);

    chain_->SetBranchStatus("L2Muon_trkPt", 1);
    chain_->SetBranchAddress("L2Muon_trkPt", &L2Muon_trkPt);

    chain_->SetBranchStatus("nTkMuon", 1);
    chain_->SetBranchAddress("nTkMuon", &nTkMuon);

    chain_->SetBranchStatus("TkMuon_pt", 1);
    chain_->SetBranchAddress("TkMuon_pt", &TkMuon_pt);

    chain_->SetBranchStatus("TkMuon_eta", 1);
    chain_->SetBranchAddress("TkMuon_eta", &TkMuon_eta);

    chain_->SetBranchStatus("TkMuon_phi", 1);
    chain_->SetBranchAddress("TkMuon_phi", &TkMuon_phi);

    chain_->SetBranchStatus("TkMuon_charge", 1);
    chain_->SetBranchAddress("TkMuon_charge", &TkMuon_charge);

    chain_->SetBranchStatus("TkMuon_trkPt", 1);
    chain_->SetBranchAddress("TkMuon_trkPt", &TkMuon_trkPt);

    chain_->SetBranchStatus("nL1Muon", 1);
    chain_->SetBranchAddress("nL1Muon", &nL1Muon);

    chain_->SetBranchStatus("L1Muon_pt", 1);
    chain_->SetBranchAddress("L1Muon_pt", &L1Muon_pt);

    chain_->SetBranchStatus("L1Muon_eta", 1);
    chain_->SetBranchAddress("L1Muon_eta", &L1Muon_eta);

    chain_->SetBranchStatus("L1Muon_phi", 1);
    chain_->SetBranchAddress("L1Muon_phi", &L1Muon_phi);

    chain_->SetBranchStatus("L1Muon_charge", 1);
    chain_->SetBranchAddress("L1Muon_charge", &L1Muon_charge);

    chain_->SetBranchStatus("L1Muon_quality", 1);
    chain_->SetBranchAddress("L1Muon_quality", &L1Muon_quality);
  }

  void TurnOnBranches_IterL3Muon()
  {
    chain_->SetBranchStatus("nIterL3OI", 1);
    chain_->SetBranchAddress("nIterL3OI", &nIterL3OI);

    chain_->SetBranchStatus("iterL3OI_inner_pt", 1);
    chain_->SetBranchAddress("iterL3OI_inner_pt", &iterL3OI_inner_pt);

    chain_->SetBranchStatus("iterL3OI_inner_eta", 1);
    chain_->SetBranchAddress("iterL3OI_inner_eta", &iterL3OI_inner_eta);

    chain_->SetBranchStatus("iterL3OI_inner_phi", 1);
    chain_->SetBranchAddress("iterL3OI_inner_phi", &iterL3OI_inner_phi);

    chain_->SetBranchStatus("iterL3OI_inner_charge", 1);
    chain_->SetBranchAddress("iterL3OI_inner_charge", &iterL3OI_inner_charge);

    chain_->SetBranchStatus("iterL3OI_outer_pt", 1);
    chain_->SetBranchAddress("iterL3OI_outer_pt", &iterL3OI_outer_pt);

    chain_->SetBranchStatus("iterL3OI_outer_eta", 1);
    chain_->SetBranchAddress("iterL3OI_outer_eta", &iterL3OI_outer_eta);

    chain_->SetBranchStatus("iterL3OI_outer_phi", 1);
    chain_->SetBranchAddress("iterL3OI_outer_phi", &iterL3OI_outer_phi);

    chain_->SetBranchStatus("iterL3OI_outer_charge", 1);
    chain_->SetBranchAddress("iterL3OI_outer_charge", &iterL3OI_outer_charge);

    chain_->SetBranchStatus("iterL3OI_global_pt", 1);
    chain_->SetBranchAddress("iterL3OI_global_pt", &iterL3OI_global_pt);

    chain_->SetBranchStatus("iterL3OI_global_eta", 1);
    chain_->SetBranchAddress("iterL3OI_global_eta", &iterL3OI_global_eta);

    chain_->SetBranchStatus("iterL3OI_global_phi", 1);
    chain_->SetBranchAddress("iterL3OI_global_phi", &iterL3OI_global_phi);

    chain_->SetBranchStatus("iterL3OI_global_charge", 1);
    chain_->SetBranchAddress("iterL3OI_global_charge", &iterL3OI_global_charge);

    chain_->SetBranchStatus("nIterL3IOFromL2", 1);
    chain_->SetBranchAddress("nIterL3IOFromL2", &nIterL3IOFromL2);

    chain_->SetBranchStatus("iterL3IOFromL2_inner_pt", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_inner_pt", &iterL3IOFromL2_inner_pt);

    chain_->SetBranchStatus("iterL3IOFromL2_inner_eta", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_inner_eta", &iterL3IOFromL2_inner_eta);

    chain_->SetBranchStatus("iterL3IOFromL2_inner_phi", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_inner_phi", &iterL3IOFromL2_inner_phi);

    chain_->SetBranchStatus("iterL3IOFromL2_inner_charge", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_inner_charge", &iterL3IOFromL2_inner_charge);

    chain_->SetBranchStatus("iterL3IOFromL2_outer_pt", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_outer_pt", &iterL3IOFromL2_outer_pt);

    chain_->SetBranchStatus("iterL3IOFromL2_outer_eta", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_outer_eta", &iterL3IOFromL2_outer_eta);

    chain_->SetBranchStatus("iterL3IOFromL2_outer_phi", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_outer_phi", &iterL3IOFromL2_outer_phi);

    chain_->SetBranchStatus("iterL3IOFromL2_outer_charge", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_outer_charge", &iterL3IOFromL2_outer_charge);

    chain_->SetBranchStatus("iterL3IOFromL2_global_pt", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_global_pt", &iterL3IOFromL2_global_pt);

    chain_->SetBranchStatus("iterL3IOFromL2_global_eta", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_global_eta", &iterL3IOFromL2_global_eta);

    chain_->SetBranchStatus("iterL3IOFromL2_global_phi", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_global_phi", &iterL3IOFromL2_global_phi);

    chain_->SetBranchStatus("iterL3IOFromL2_global_charge", 1);
    chain_->SetBranchAddress("iterL3IOFromL2_global_charge", &iterL3IOFromL2_global_charge);

    chain_->SetBranchStatus("nIterL3FromL2", 1);
    chain_->SetBranchAddress("nIterL3FromL2", &nIterL3FromL2);

    chain_->SetBranchStatus("iterL3FromL2_inner_pt", 1);
    chain_->SetBranchAddress("iterL3FromL2_inner_pt", &iterL3FromL2_inner_pt);

    chain_->SetBranchStatus("iterL3FromL2_inner_eta", 1);
    chain_->SetBranchAddress("iterL3FromL2_inner_eta", &iterL3FromL2_inner_eta);

    chain_->SetBranchStatus("iterL3FromL2_inner_phi", 1);
    chain_->SetBranchAddress("iterL3FromL2_inner_phi", &iterL3FromL2_inner_phi);

    chain_->SetBranchStatus("iterL3FromL2_inner_charge", 1);
    chain_->SetBranchAddress("iterL3FromL2_inner_charge", &iterL3FromL2_inner_charge);

    chain_->SetBranchStatus("iterL3FromL2_outer_pt", 1);
    chain_->SetBranchAddress("iterL3FromL2_outer_pt", &iterL3FromL2_outer_pt);

    chain_->SetBranchStatus("iterL3FromL2_outer_eta", 1);
    chain_->SetBranchAddress("iterL3FromL2_outer_eta", &iterL3FromL2_outer_eta);

    chain_->SetBranchStatus("iterL3FromL2_outer_phi", 1);
    chain_->SetBranchAddress("iterL3FromL2_outer_phi", &iterL3FromL2_outer_phi);

    chain_->SetBranchStatus("iterL3FromL2_outer_charge", 1);
    chain_->SetBranchAddress("iterL3FromL2_outer_charge", &iterL3FromL2_outer_charge);

    chain_->SetBranchStatus("iterL3FromL2_global_pt", 1);
    chain_->SetBranchAddress("iterL3FromL2_global_pt", &iterL3FromL2_global_pt);

    chain_->SetBranchStatus("iterL3FromL2_global_eta", 1);
    chain_->SetBranchAddress("iterL3FromL2_global_eta", &iterL3FromL2_global_eta);

    chain_->SetBranchStatus("iterL3FromL2_global_phi", 1);
    chain_->SetBranchAddress("iterL3FromL2_global_phi", &iterL3FromL2_global_phi);

    chain_->SetBranchStatus("iterL3FromL2_global_charge", 1);
    chain_->SetBranchAddress("iterL3FromL2_global_charge", &iterL3FromL2_global_charge);

    chain_->SetBranchStatus("nIterL3IOFromL1", 1);
    chain_->SetBranchAddress("nIterL3IOFromL1", &nIterL3IOFromL1);

    chain_->SetBranchStatus("iterL3IOFromL1_pt", 1);
    chain_->SetBranchAddress("iterL3IOFromL1_pt", &iterL3IOFromL1_pt);

    chain_->SetBranchStatus("iterL3IOFromL1_eta", 1);
    chain_->SetBranchAddress("iterL3IOFromL1_eta", &iterL3IOFromL1_eta);

    chain_->SetBranchStatus("iterL3IOFromL1_phi", 1);
    chain_->SetBranchAddress("iterL3IOFromL1_phi", &iterL3IOFromL1_phi);

    chain_->SetBranchStatus("iterL3IOFromL1_charge", 1);
    chain_->SetBranchAddress("iterL3IOFromL1_charge", &iterL3IOFromL1_charge);

    chain_->SetBranchStatus("nIterL3MuonNoID", 1);
    chain_->SetBranchAddress("nIterL3MuonNoID", &nIterL3MuonNoID);

    chain_->SetBranchStatus("iterL3MuonNoID_pt", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_pt", &iterL3MuonNoID_pt);

    chain_->SetBranchStatus("iterL3MuonNoID_eta", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_eta", &iterL3MuonNoID_eta);

    chain_->SetBranchStatus("iterL3MuonNoID_phi", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_phi", &iterL3MuonNoID_phi);

    chain_->SetBranchStatus("iterL3MuonNoID_charge", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_charge", &iterL3MuonNoID_charge);

    chain_->SetBranchStatus("iterL3MuonNoID_isGLB", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_isGLB", &iterL3MuonNoID_isGLB);

    chain_->SetBranchStatus("iterL3MuonNoID_isSTA", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_isSTA", &iterL3MuonNoID_isSTA);

    chain_->SetBranchStatus("iterL3MuonNoID_isTRK", 1);
    chain_->SetBranchAddress("iterL3MuonNoID_isTRK", &iterL3MuonNoID_isTRK);    
  }
};

};