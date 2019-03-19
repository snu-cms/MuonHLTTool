// -- Example file on how to read the object information in the ntuples
// -- command: root -l -b -q PrintObject.cxx

#include <Include/MuonHLTTool.h>

class Example
{
public:
  vector<TString> vec_ntuplePath_;

  Example()
  {

  }

  void AddNtuplePath(TString ntuplePath)
  {
    vec_ntuplePath_.push_back( ntuplePath );
  }

  void Run()
  {
    StartTimer();

    TChain *chain = new TChain("ntupler/ntuple");
    for(const auto& ntuplePath : vec_ntuplePath_ )
      chain->Add( ntuplePath );

    MuonHLT::NtupleHandle *ntuple = new MuonHLT::NtupleHandle(chain);
    // -- turn-on branches if you need - e.g. TurnOnBranches_Muon(), TurnOnBranches_GenParticle()
    ntuple->TurnOnBranches_HLTMuon();

    Int_t nEvent = chain->GetEntries();
    for(Int_t i_ev=0; i_ev<nEvent; i_ev++)
    {
      chain->GetEvent(i_ev);
      cout << "*** " << i_ev << " th event ***" << endl;

      cout << "[List of fired trigger]" << endl;
      for(const auto& firedTrigger: *ntuple->vec_firedTrigger )
      {
        cout << "  " << firedTrigger << endl;
      }

      cout << "[List of my-fired trigger (from HLT rerunning)]" << endl;
      for(const auto& firedTrigger: *ntuple->vec_myFiredTrigger )
      {
        cout << "  " << firedTrigger << endl;
      }

      vector<MuonHLT::HLTObject> vec_HLTObject = MuonHLT::GetAllHLTObject(ntuple);
      cout << "[List of all HLT objects]" << endl;
      for(const auto& HLTObject : vec_HLTObject )
      {
        cout << "  (filterName, pt, eta, phi) = ("
        << HLTObject.filterName << ", "
        << HLTObject.pt         << ", "
        << HLTObject.eta        << ", "
        << HLTObject.phi        << ")" << endl;
      }

      vector<MuonHLT::MYHLTObject> vec_MYHLTObject = MuonHLT::GetAllMYHLTObject(ntuple);
      cout << "[List of all MYHLT objects (from HLT rerunning)]" << endl;
      for(const auto& MYHLTObject : vec_MYHLTObject )
      {
        cout << "  (filterName, pt, eta, phi) = ("
        << MYHLTObject.filterName << ", "
        << MYHLTObject.pt         << ", "
        << MYHLTObject.eta        << ", "
        << MYHLTObject.phi        << ")" << endl;
      }

      vector<MuonHLT::L1Muon> vec_L1Muon = MuonHLT::GetAllL1Muon(ntuple);
      cout << "[List of all L1 muons]" << endl;
      for(const auto& l1mu : vec_L1Muon )
      {
        cout << "  (pt, eta, phi, charge, quality) = ("
        << l1mu.pt      << ", "
        << l1mu.eta     << ", "
        << l1mu.phi     << ", "
        << l1mu.charge  << ", "
        << l1mu.quality << ")" << endl;
      }

      vector<MuonHLT::L2Muon> vec_L2Muon = MuonHLT::GetAllL2Muon(ntuple);
      cout << "[List of all L2 muons]" << endl;
      for(const auto& l2mu : vec_L2Muon )
      {
        cout << "  (pt, eta, phi, charge, TrkPt) = ("
        << l2mu.pt      << ", "
        << l2mu.eta     << ", "
        << l2mu.phi     << ", "
        << l2mu.charge  << ", "
        << l2mu.trkPt   << ")" << endl;
      }

      cout << endl;

      if( i_ev > 100 ) break;
    }

    PrintRunTime();
  }

private:
  TStopwatch timer_;

  void StartTimer()
  {
    timer_.Start();
  }

  void PrintRunTime()
  {
    Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << "************************************************" << endl;
    cout << "Total real time: " << realTime << " (seconds)" << endl;
    cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
    cout << "  CPU time / real time = " << cpuTime / realTime << endl;
    cout << "************************************************" << endl;
  }
};

void PrintObject()
{
  TString analyzerPath = gSystem->Getenv("MUONHLT_ANALYZER_PATH");
  TString ntuplePath = analyzerPath+"/Example/ExampleNtuple_ZMuMu_M50to120.root";

  Example* example = new Example();
  example->AddNtuplePath(ntuplePath);
  example->Run();  
}