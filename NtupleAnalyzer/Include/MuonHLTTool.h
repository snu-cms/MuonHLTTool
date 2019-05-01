#include <Include/Object.h>
#include <iostream>

namespace MuonHLT
{
// -- some useful variables, functions, classes ...

vector<MuonHLT::HLTObject> GetAllHLTObject(NtupleHandle *ntuple, TString filterName = "", Double_t minPt = -1)
{
  vector<MuonHLT::HLTObject> vec_HLTObj;
  Int_t nHLTObj = (Int_t)ntuple->vec_filterName->size();
  for(Int_t i_obj=0; i_obj<nHLTObj; i_obj++)
  {
    MuonHLT::HLTObject hltObject(ntuple, i_obj);

    if( hltObject.pt < minPt ) continue;
    if( filterName != "" && hltObject.filterName != filterName ) continue;

    // -- if filterName is not specified: save all HLT object
    // -- if filterName is speficied: save HLT object with same filter name only
    vec_HLTObj.push_back( hltObject );
  }

  return vec_HLTObj;
}

vector<MuonHLT::MYHLTObject> GetAllMYHLTObject(NtupleHandle *ntuple, TString filterName = "", Double_t minPt = -1)
{
  vector<MuonHLT::MYHLTObject> vec_MYHLTObj;
  Int_t nHLTObj = (Int_t)ntuple->vec_myFilterName->size();
  for(Int_t i_obj=0; i_obj<nHLTObj; i_obj++)
  {
    MuonHLT::MYHLTObject myHLTObject(ntuple, i_obj);

    if( myHLTObject.pt < minPt ) continue;
    if( filterName != "" && myHLTObject.filterName != filterName ) continue;

    // -- if filterName is not specified: save all HLT object
    // -- if filterName is speficied: save HLT object with same filter name only
    vec_MYHLTObj.push_back( myHLTObject );
  }

  return vec_MYHLTObj;
}

vector<MuonHLT::L1Muon> GetAllL1Muon(NtupleHandle* ntuple, Double_t minPt = -1 )
{
  vector<MuonHLT::L1Muon> vec_L1Muon;
  for(Int_t i_L1=0; i_L1<ntuple->nL1Muon; i_L1++)
  {
    MuonHLT::L1Muon l1mu(ntuple, i_L1);
    if( l1mu.pt > minPt )
      vec_L1Muon.push_back( l1mu );
  }

  return vec_L1Muon;
}

vector<MuonHLT::L2Muon> GetAllL2Muon(NtupleHandle* ntuple, Double_t minPt = -1 )
{
  vector<MuonHLT::L2Muon> vec_muon;
  for(Int_t i_obj=0; i_obj<ntuple->nL2Muon; i_obj++)
  {
    MuonHLT::L2Muon mu(ntuple, i_obj);
    if( mu.pt > minPt )
      vec_muon.push_back( mu );
  }

  return vec_muon;
}

Bool_t dRMatching( TLorentzVector vecP_ref, vector<TLorentzVector> vec_vecP, Double_t minDR )
{
  bool flag = kFALSE;

  Int_t nObj = (Int_t)vec_vecP.size();
  for(const auto& vecP_target: vec_vecP )
  {
    Double_t dR = vecP_ref.DeltaR( vecP_target );
    if( dR < minDR )
    {
      flag = kTRUE;
      break;
    }
  }

  return flag;
}

Bool_t dRMatching_L1Muon( TLorentzVector vecP_ref, NtupleHandle* ntuple, Double_t minPt, Double_t minDR )
{
  vector<MuonHLT::L1Muon> vec_L1Muon = MuonHLT::GetAllL1Muon( ntuple, minPt );
  vector<TLorentzVector> vec_vecP_L1Muon;
  for(const auto& L1Muon : vec_L1Muon )
    vec_vecP_L1Muon.push_back( L1Muon.vecP );

  return MuonHLT::dRMatching( vecP_ref, vec_vecP_L1Muon, minDR );
}

Bool_t dRMatching_HLTObj( TLorentzVector vecP_ref, MuonHLT::NtupleHandle* ntuple, TString filterName, Double_t minDR )
{
  vector<MuonHLT::HLTObject> vec_HLTObj = MuonHLT::GetAllHLTObject(ntuple, filterName);
  vector<TLorentzVector> vec_vecP_HLTObj;
  for(const auto& HLTObj : vec_HLTObj )
    vec_vecP_HLTObj.push_back( HLTObj.vecP );

  return MuonHLT::dRMatching( vecP_ref, vec_vecP_HLTObj, minDR );
}

// -- same with above function, but take MYHLT object
Bool_t dRMatching_MYHLTObj( TLorentzVector vecP_ref, MuonHLT::NtupleHandle* ntuple, TString filterName, Double_t minDR )
{
  vector<MuonHLT::MYHLTObject> vec_HLTObj = MuonHLT::GetAllMYHLTObject(ntuple, filterName);
  vector<TLorentzVector> vec_vecP_HLTObj;
  for(const auto& HLTObj : vec_HLTObj )
    vec_vecP_HLTObj.push_back( HLTObj.vecP );

  return MuonHLT::dRMatching( vecP_ref, vec_vecP_HLTObj, minDR );
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
      cout << endl;

    if ( x % (n/r +1) != 0 ) return;

 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";
 
    for (int x=c; x<w; x++) cout << " ";
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
  cout << "]\r" << flush;

}


}; // -- end of namespace
