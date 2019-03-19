#pragma once

#include <TH1D.h>
#include <TH2D.h>
#include <TColor.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TFile.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <THStack.h>

#include <vector>

namespace PlotTool {

TH1D* Get_Hist(TString fileName, TString histName, TString histName_new = "" )
{
  TH1::AddDirectory(kFALSE);

  TFile* f_input = TFile::Open( fileName );
  TH1D* h_temp = (TH1D*)f_input->Get(histName)->Clone();
  if( histName_new != "" )
    h_temp->SetName( histName_new );

  f_input->Close();

  return h_temp;
}

TH2D* Get_Hist2D(TString fileName, TString histName, TString histName_new = "" )
{
  TH1::AddDirectory(kFALSE);

  TFile* f_input = TFile::Open( fileName );
  TH2D* h_temp = (TH2D*)f_input->Get(histName)->Clone();
  if( histName_new != "" )
    h_temp->SetName( histName_new );

  f_input->Close();

  return h_temp;
}


TGraphAsymmErrors* Get_Graph(TString fileName, TString graphName, TString graphName_New = "" )
{
  TFile *f_input = TFile::Open( fileName );
  TGraphAsymmErrors* g_temp = (TGraphAsymmErrors*)f_input->Get(graphName)->Clone();
  if( graphName_New != "" )
    g_temp->SetName( graphName_New );

  f_input->Close();

  return g_temp;
}

void SetLegend( TLegend *& legend, Double_t xMin = 0.75, Double_t yMin = 0.75, Double_t xMax = 0.95, Double_t yMax = 0.95 )
{
  legend = new TLegend( xMin, yMin, xMax, yMax );
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont( 62 );
}

void SetAxis_SinglePad( TAxis *axisX, TAxis *axisY, TString titleX, TString titleY )
{
  axisX->SetTitle( titleX );
  axisX->SetTitleFont(42);
  axisX->SetTitleSize(0.05);
  axisX->SetTitleOffset(1.1);

  axisX->SetLabelFont(42);
  axisX->SetLabelSize(0.04);
  axisX->SetNoExponent();
  axisX->SetMoreLogLabels();

  axisY->SetTitle( titleY );
  axisY->SetTitleFont(42);
  axisY->SetTitleSize(0.05);
  axisY->SetTitleOffset(1.2);

  axisY->SetLabelFont(42);
  axisY->SetLabelSize(0.04);
}

void SetAxis_TopPad( TAxis *axisX, TAxis *axisY, TString titleY )
{
  axisX->SetLabelFont(42);
  axisX->SetLabelSize(0.000);
  axisX->SetTitleSize(0.000);

  axisY->SetTitle( titleY );
  axisY->SetTitleFont(42);
  axisY->SetTitleSize(0.05);
  axisY->SetTitleOffset(1.25);

  axisY->SetLabelFont(42);
  axisY->SetLabelSize(0.04);
}

void SetAxis_BottomPad( TAxis *axisX, TAxis *axisY, TString titleX, TString titleY)
{

  axisX->SetTitle( titleX );
  axisX->SetTitleFont(42);
  axisX->SetTitleSize( 0.2 );
  axisX->SetTitleOffset( 0.85 );

  axisX->SetLabelFont(42);
  axisX->SetLabelSize(0.13);
  axisX->SetLabelOffset(0.01);
  axisX->SetLabelColor(1);
  axisX->SetMoreLogLabels();
  axisX->SetNoExponent();


  axisY->SetTitle( titleY );
  axisY->SetTitleFont(42);
  axisY->SetTitleSize(0.12);
  axisY->SetTitleOffset( 0.55 );

  axisY->SetLabelFont(42);
  axisY->SetLabelSize( 0.10 );
}

void DrawLine( TF1*& f_line, Int_t color = kRed )
{
  f_line = new TF1("f_line", "1", -10000, 10000);
  f_line->SetLineColor(color);
  f_line->SetLineWidth(1);
  f_line->Draw("PSAME");
}

TH1D* DivideEachBin_ByBinWidth( TH1D* h, TString HistName = "" )
{
  TH1D* h_return = (TH1D*)h->Clone();
  if( HistName != "" )
    h_return->SetName(HistName);

  Int_t nBin = h->GetNbinsX();
  for(Int_t i=0; i<nBin; i++)
  {
    Int_t i_bin = i+1;
    Double_t Entry_before = h->GetBinContent(i_bin);
    Double_t Error_before = h->GetBinError(i_bin);
    Double_t BinWidth = h->GetBinWidth(i_bin);

    Double_t Entry_after = Entry_before / BinWidth;
    Double_t Error_after = Error_before / BinWidth;

    h_return->SetBinContent(i_bin, Entry_after);
    h_return->SetBinError(i_bin, Error_after);
  }

  return h_return;
}

TH1D* MultiplyEachBin_byBinWidth( TH1D* h, TString HistName = "" )
{
  TH1D* h_return = (TH1D*)h->Clone();
  if( HistName != "" )
    h_return->SetName(HistName);

  Int_t nBin = h->GetNbinsX();
  for(Int_t i=0; i<nBin; i++)
  {
    Int_t i_bin = i+1;
    Double_t Entry_before = h->GetBinContent(i_bin);
    Double_t Error_before = h->GetBinError(i_bin);
    Double_t BinWidth = h->GetBinWidth(i_bin);

    Double_t Entry_after = Entry_before * BinWidth;
    Double_t Error_after = Error_before * BinWidth;

    h_return->SetBinContent(i_bin, Entry_after);
    h_return->SetBinError(i_bin, Error_after);
  }

  return h_return;
}

Bool_t IsRatio1( TH1D* h1, TH1D* h2)
{
  Bool_t isRatio1 = kTRUE;

  TString h1Name = h1->GetName();
  TString h2Name = h2->GetName();
  printf("[IsRatio1] Check %s / %s == 1 for all bin\n", h1Name.Data(), h2Name.Data());

  Int_t nBin1 = h1->GetNbinsX();
  Int_t nBin2 = h2->GetNbinsX();
  if( nBin1 != nBin2 )
  {
    printf("(nBin1 = %d, nBin2 = %d) is different! ... need to check\n", nBin1, nBin2);
    return kFALSE;
  }

  for(Int_t i=0; i<nBin1; i++)
  {
    Int_t i_bin = i+1;

    Double_t content1 = h1->GetBinContent(i_bin);
    Double_t content2 = h2->GetBinContent(i_bin);
    Double_t ratio = content1 / content2;

    if( fabs(ratio - 1) > 1e-5 )
    {
      printf("[%02d bin is deviated from 1]\n", i_bin);
      printf("   Bin content1 = %lf\n", content1);
      printf("   Bin content2 = %lf\n", content2);
      printf("          Ratio = %lf\n\n", ratio);
      isRatio1 = kFALSE;
    }
  }

  TString isSuccess = isRatio1 ? "success" : "fail";

  printf("\n");
  printf("[IsRatio1] Checking is finished\n");
  printf("   Result: %s\n\n", isSuccess.Data());

  return isRatio1;
}

TH1D* Convert_GraphToHist( TGraphAsymmErrors *g )
{
  const Int_t nBin = g->GetN();
  Double_t *BinEdges = new Double_t[nBin+1];
  Double_t *value = new Double_t[nBin];
  Double_t *error = new Double_t[nBin];

  for(Int_t i=0; i<nBin; i++)
  {
    Double_t x, y;
    g->GetPoint(i, x, y);

    // -- make BinEdges array -- //
    Double_t ErrX_Low = g->GetErrorXlow(i);
    Double_t ErrX_High = g->GetErrorXhigh(i);

    if( i == nBin-1 )
    {
      BinEdges[i] = x - ErrX_Low;
      BinEdges[i+1] = x + ErrX_High;
    }
    else
      BinEdges[i] = x - ErrX_Low;


    // -- store graph information -- //
    value[i] = y;

    Double_t ErrY_Low = g->GetErrorYlow(i);
    Double_t ErrY_High = g->GetErrorYhigh(i);

    // -- take the larger one -- //
    error[i] = ErrY_Low > ErrY_High ? ErrY_Low : ErrY_High;
  }

  TString GraphName = g->GetName();
  TH1D* h_temp = new TH1D( "h_"+GraphName, "", nBin, BinEdges );

  // -- fill this histogram using graph information -- //
  for(Int_t i=0; i<nBin; i++)
  {
    Int_t i_bin = i+1;
    h_temp->SetBinContent( i_bin, value[i] );
    h_temp->SetBinError( i_bin, error[i] );
  }

  return h_temp;
}

struct HistInfo
{
  TH1D* h;
  TString legend;
  Int_t color;
};

struct GraphInfo
{
  TGraphAsymmErrors* g;
  TString legend;
  Int_t color;
};

struct LatexInfo
{
  Double_t x;
  Double_t y;
  TString content;
};

class CanvasBase
{
public:
  TCanvas* c_;
  TString canvasName_;

  TString titleX_;
  TString titleY_;

  Bool_t isLogX_;
  Bool_t isLogY_;

  Double_t legendMinX_;
  Double_t legendMaxX_;
  Double_t legendMinY_;
  Double_t legendMaxY_;

  Bool_t setLegendColumn_;
  Int_t nLegendColumn_;

  Double_t minX_;
  Double_t maxX_;
  Bool_t setRangeX_;

  Double_t minY_;
  Double_t maxY_;
  Bool_t setRangeY_;

  // -- latex (CMS Preliminary, lumi. info, etc.)
  TLatex latex_;
  Bool_t setLatexCMSPre_;
  Bool_t setLatexLumiEnergy_;
  Double_t lumi_;
  Int_t energy_;
  Bool_t setLatexCMSSim_;

  // -- additional latex info.
  vector<LatexInfo> latexInfos_;
  Bool_t setLatexInfo_;

  // -- for the canvas with ratio plot
  TPad* topPad_;
  TPad* bottomPad_;

  TString titleRatio_;

  Double_t minRatio_;
  Double_t maxRatio_;
  Bool_t setRangeRatio_;

  // -- for auto-adjustment of y-axis
  Bool_t setAutoRangeY_;

  Bool_t doRemoveZeroPoint_;

  CanvasBase()
  {
    Init();
  }

  CanvasBase(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE )
  {
    canvasName_ = canvasName_;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void SetTitle( TString titleX, TString titleY )
  {
    titleX_ = titleX;
    titleY_ = titleY;
  }

  void SetTitle( TString titleX, TString titleY, TString titleRatio )
  {
    titleX_ = titleX;
    titleY_ = titleY;
    titleRatio_ = titleRatio;
  }

  void SetLegendPosition( Double_t minX, Double_t minY, Double_t maxX, Double_t maxY )
  {
    legendMinX_ = minX;
    legendMinY_ = minY;
    legendMaxX_ = maxX;
    legendMaxY_ = maxY;
  }

  void SetLegendColumn( Int_t nColumn )
  {
    setLegendColumn_ = kTRUE;
    nLegendColumn_ = nColumn;
  }

  void SetRangeX( Double_t min, Double_t max )
  {
    minX_ = min;
    maxX_ = max;
    setRangeX_ = kTRUE;
  }

  void SetRangeY( Double_t min, Double_t max )
  {
    minY_ = min;
    maxY_ = max;
    setRangeY_ = kTRUE;
  }

  void SetRangeRatio( Double_t min, Double_t max )
  {
    minRatio_ = min;
    maxRatio_ = max;
    setRangeRatio_ = kTRUE;
  }

  void Latex_CMSPre()
  {
    setLatexCMSPre_ = kTRUE;
  }

  void Latex_CMSPre(Double_t lumi, Int_t energy)
  {
    Latex_CMSPre();
    setLatexLumiEnergy_ = kTRUE;
    lumi_ = lumi;
    energy_ = energy;
  }

  void Latex_CMSSim()
  {
    setLatexCMSSim_ = kTRUE;
  }

  void RegisterLatex( Double_t x, Double_t y, TString content )
  {
    setLatexInfo_ = kTRUE;
    LatexInfo latexInfo{x, y, content};
    latexInfos_.push_back( latexInfo );
  }

  // -- implemented later
  virtual void Draw( TString drawOp )
  {

  }

  void Init()
  {
    canvasName_ = "undefined";
    isLogX_ = kFALSE;
    isLogY_ = kFALSE;

    titleX_ = "undefined";
    titleY_ = "undefined";

    legendMinX_ = 0.50;
    legendMinY_ = 0.70;
    legendMaxX_ = 0.95;
    legendMaxY_ = 0.95;

    setLegendColumn_ = kFALSE;
    nLegendColumn_ = 1;

    setRangeX_ = kFALSE;
    minX_ = 0;
    maxX_ = 0;

    setRangeY_ = kFALSE;
    minY_ = 0;
    maxY_ = 0;

    setLatexCMSPre_ = kFALSE;
    setLatexLumiEnergy_ = kFALSE;
    lumi_ = -999;
    energy_ = -999;
    setLatexCMSSim_ = kFALSE;
    setLatexInfo_ = kFALSE;

    // -- for the canvas with ratio plot
    topPad_ = NULL;
    bottomPad_ = NULL;

    titleRatio_ = "undefined";

    setRangeRatio_ = kFALSE;
    minRatio_ = 0;
    maxRatio_ = 2.5;

    setAutoRangeY_ = kFALSE;

    doRemoveZeroPoint_ = kFALSE;
  }

  void DrawLatex_CMSPre()
  {
    latex_.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Preliminary}}}");
  }

  void DrawLatex_CMSPreLumiEnergy()
  {
    DrawLatex_CMSPre();
    latex_.DrawLatexNDC(0.70, 0.96, "#font[42]{#scale[0.7]{"+TString::Format("%.1lf fb^{-1} (%d TeV)", lumi_, energy_)+"}}");
  }

  void DrawLatex_CMSSim()
  {
    latex_.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Simulation}}}");
    latex_.DrawLatexNDC(0.85, 0.96, "#font[42]{#scale[0.7]{13 TeV}}");
  }

  void SetCanvas_Square()
  {
    c_ = new TCanvas(canvasName_, "", 800, 800);
    c_->cd();
    
    c_->SetTopMargin(0.05);
    c_->SetLeftMargin(0.13);
    c_->SetRightMargin(0.045);
    c_->SetBottomMargin(0.13);

    if( isLogX_ )
      c_->SetLogx();
    if( isLogY_ )
      c_->SetLogy();
  }

  void SetCanvas_Ratio()
  {
    c_ = new TCanvas(canvasName_, "", 800, 800);
    c_->cd();

    topPad_ = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    topPad_->Draw();
    topPad_->cd();

    topPad_->SetTopMargin(0.05);
    topPad_->SetLeftMargin(0.13);
    topPad_->SetRightMargin(0.045);
    topPad_->SetBottomMargin(0.3);

    if( isLogX_ ) topPad_->SetLogx();
    if( isLogY_ ) topPad_->SetLogy();

    c_->cd();
    bottomPad_ = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.29 );
    bottomPad_->Draw();
    bottomPad_->cd();
    bottomPad_->SetGridx();
    bottomPad_->SetGridy();
    bottomPad_->SetTopMargin(0.05);
    bottomPad_->SetBottomMargin(0.4);
    bottomPad_->SetRightMargin(0.045);
    bottomPad_->SetLeftMargin(0.13);

    if( isLogX_ ) bottomPad_->SetLogx();
  }

  void DrawLatexAll()
  {
    if( setLatexCMSPre_ )
    {
      if( setLatexLumiEnergy_ ) DrawLatex_CMSPreLumiEnergy();
      else                      DrawLatex_CMSPre();
    }

    if( setLatexCMSSim_ ) DrawLatex_CMSSim();

    if( setLatexInfo_ )
    {
      for( auto latexInfo : latexInfos_ )
        latex_.DrawLatexNDC( latexInfo.x, latexInfo.y, latexInfo.content );
    }
  }

  // -- for auto adjustment of Y-range
  void SetAutoRangeY(vector<TH1D*> vec_hist)
  {
    Double_t globalMin = 9999;
    Double_t globalMax = -9999;
    for(const auto& h : vec_hist )
    {
      Double_t localMin = h->GetMinimum();
      Double_t localMax = h->GetMaximum();
      if( localMin < globalMin ) globalMin = localMin;
      if( localMax > globalMax ) globalMax = localMax;
    }

    minY_ = globalMin > 0 ? 0 : globalMin * 1.3;
    maxY_ = globalMax * 1.3;

    // -- TO-DO: multiplication factor adjustment according to isLogY
    if( minY_ == 0 && isLogY_ ) minY_ = 0.5;
    if( isLogY_ ) maxY_ = globalMax * 1e2;
  }

  void RemoveZeroPoint(Bool_t value = kTRUE )
  {
    doRemoveZeroPoint_ = value;
  }

  void RemoveZeroPointInGraph(TGraphAsymmErrors* g)
  {
    for(Int_t i=0; i<g->GetN(); i++)
    {
      Double_t x, y;
      g->GetPoint(i, x, y);
      if( y == 0 )
      {
        g->RemovePoint(i);
        i = -1; // -- go back to 0 as the point is reordered after removing a point
      }
    }
  }

}; // class CanvasBase

class HistCanvas : public CanvasBase
{
public:
  vector<HistInfo> histInfos_;

  Double_t nRebin_ = 0;
  Bool_t setRebin_ = kFALSE;

  HistCanvas()
  {
    // -- member variables are initialized by Init() in CanvasBase()
  }

  HistCanvas(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE ): HistCanvas()
  {
    canvasName_ = canvasName;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void Register( TH1D* h, TString legend, Int_t color  )
  {
    HistInfo histInfo{ (TH1D*)h->Clone(), legend, color };
    histInfos_.push_back( histInfo );
  }

  void SetRebin( Int_t n )
  {
    nRebin_ = n;
    setRebin_ = kTRUE;
  }

  void Draw( TString drawOp = "EPSAME" )
  {
    if( !drawOp.Contains("SAME") ) drawOp = drawOp + "SAME";

    TLegend *legend;
    PlotTool::SetLegend( legend, legendMinX_, legendMinY_, legendMaxX_, legendMaxY_ );
    if( setLegendColumn_ ) legend->SetNColumns(nLegendColumn_);

    // -- draw canvas
    SetCanvas_Square();

    c_->cd();

    Int_t nHist = histInfos_.size();
    for(Int_t i=0; i<nHist; i++)
    {
      TH1D*& h = histInfos_[i].h;
      TString legendName = histInfos_[i].legend;
      Int_t color = histInfos_[i].color;

      if( setRebin_ ) h->Rebin( nRebin_ );

      h->Draw(drawOp);
      h->SetStats(kFALSE);
      h->SetMarkerStyle(20);
      h->SetMarkerColor(color);
      h->SetLineColor(color);
      h->SetFillColorAlpha(kWhite, 0); 
      h->SetTitle("");

      if( i == 0 ) PlotTool::SetAxis_SinglePad( h->GetXaxis(), h->GetYaxis(), titleX_, titleY_ );
      if( setRangeX_ ) h->GetXaxis()->SetRangeUser( minX_, maxX_ );
      if( setRangeY_ ) h->GetYaxis()->SetRangeUser( minY_, maxY_ );

      legend->AddEntry( h, legendName );
    }

    legend->Draw();

    DrawLatexAll();

    c_->SaveAs(".pdf");
  }
}; // -- class HistCanvas

class HistCanvaswRatio: public HistCanvas
{
public:
  vector<HistInfo> histInfoRatios_;


  HistCanvaswRatio()
  {
    // -- member variables are initialized by Init() in CanvasBase()
  }

  HistCanvaswRatio(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE ): HistCanvaswRatio()
  {
    canvasName_ = canvasName;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void Draw( TString drawOp = "EPSAME" )
  {
    if( !drawOp.Contains("SAME") ) drawOp = drawOp + "SAME";

    TLegend *legend;
    PlotTool::SetLegend( legend, legendMinX_, legendMinY_, legendMaxX_, legendMaxY_ );
    if( setLegendColumn_ ) legend->SetNColumns(nLegendColumn_);

    // -- draw canvas
    SetCanvas_Ratio();

    c_->cd();
    topPad_->cd();

    Int_t nHist = histInfos_.size();
    for(Int_t i=0; i<nHist; i++)
    {
      TH1D*& h = histInfos_[i].h;
      TString legendName = histInfos_[i].legend;
      Int_t color = histInfos_[i].color;

      if( setRebin_ ) h->Rebin( nRebin_ );

      h->Draw(drawOp);
      h->SetStats(kFALSE);
      h->SetMarkerStyle(20);
      h->SetMarkerColor(color);
      h->SetLineColor(color);
      h->SetFillColorAlpha(kWhite, 0); 
      h->SetTitle("");

      if( i == 0 ) PlotTool::SetAxis_TopPad( h->GetXaxis(), h->GetYaxis(), titleY_ );
      if( setRangeX_ ) h->GetXaxis()->SetRangeUser( minX_, maxX_ );
      if( setRangeY_ ) h->GetYaxis()->SetRangeUser( minY_, maxY_ );

      legend->AddEntry( h, legendName );
    }

    legend->Draw();

    DrawLatexAll();

    // -- bottom pad
    c_->cd();
    bottomPad_->cd();

    CalcRatioHist();

    Int_t nHistRatio = histInfoRatios_.size();
    for(Int_t i=0; i<nHistRatio; i++)
    {
      TH1D*& h_ratio = histInfoRatios_[i].h;
      Int_t  color   = histInfoRatios_[i].color;

      h_ratio->Draw(drawOp);
      h_ratio->SetStats(kFALSE);
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetMarkerColor(color);
      h_ratio->SetLineColor(color);
      h_ratio->SetFillColorAlpha(kWhite, 0); 
      h_ratio->SetTitle("");
      if( i == 0 ) SetAxis_BottomPad(h_ratio->GetXaxis(), h_ratio->GetYaxis(), titleX_, titleRatio_);
      if( setRangeX_ )     h_ratio->GetXaxis()->SetRangeUser( minX_, maxX_ );
      if( setRangeRatio_ ) h_ratio->GetYaxis()->SetRangeUser( minRatio_, maxRatio_ );
    }

    TF1 *f_line;
    PlotTool::DrawLine(f_line);

    c_->SaveAs(".pdf");
  }

  void CalcRatioHist()
  {
    TH1D* h_ref = histInfos_[0].h;
    h_ref->Sumw2();

    Int_t nHist = histInfos_.size();
    for(Int_t i=1; i<nHist; i++) // -- starts with 1 -- //
    {
      TH1D* h_target = (TH1D*)histInfos_[i].h->Clone();
      h_target->Sumw2();
      
      TString legend = histInfos_[i].legend;
      Int_t color = histInfos_[i].color;

      TH1D* h_ratioTemp = (TH1D*)h_ref->Clone();
      h_ratioTemp->Divide( h_target, h_ref );

      HistInfo histInfoRatio{ h_ratioTemp, legend, color };
      histInfoRatios_.push_back( histInfoRatio );
    }
  }
};

class HistStackCanvaswRatio: public HistCanvas
{
public:
  HistInfo histInfo_data_;

  THStack *hs;
  TH1D* h_ratio_dataToStack_;


  HistStackCanvaswRatio()
  {
    // -- member variables are initialized by Init() in CanvasBase()
  }

  HistStackCanvaswRatio(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE ): HistStackCanvaswRatio()
  {
    canvasName_ = canvasName;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void RegisterData( TH1D* h, TString legend, Int_t color  )
  {
    histInfo_data_.h = (TH1D*)h->Clone();
    histInfo_data_.legend = legend;
    histInfo_data_.color = color;
  }

  void Draw( TString drawOp = "EPSAME" )
  {
    if( !drawOp.Contains("SAME") ) drawOp = drawOp + "SAME";

    // -- make legend
    TLegend *legend;
    PlotTool::SetLegend( legend, legendMinX_, legendMinY_, legendMaxX_, legendMaxY_ );
    if( setLegendColumn_ ) legend->SetNColumns(nLegendColumn_);

    // -- setup data, MC stacks + add in the legend
    SetDataHistogram(legend);
    SetMCStack(legend);

    // -- initialize the canvas
    SetCanvas_Ratio();

    c_->cd();
    topPad_->cd();

    TH1D* h_format = (TH1D*)histInfo_data_.h->Clone();
    h_format->Reset("ICES");
    h_format->Draw("");
    PlotTool::SetAxis_TopPad( h_format->GetXaxis(), h_format->GetYaxis(), titleY_ );
    if( setRangeX_ ) h_format->GetXaxis()->SetRangeUser( minX_, maxX_ );

    // -- automatic y-axis range
    vector<TH1D*> vec_histAll;
    vec_histAll.push_back( histInfo_data_.h );
    for(const auto& histInfo: histInfos_ ) vec_histAll.push_back( histInfo.h );
    SetAutoRangeY(vec_histAll);
    h_format->GetYaxis()->SetRangeUser(minY_, maxY_);

    // -- but if a specific y-range is provided, then it should override
    if( setRangeY_ ) h_format->GetYaxis()->SetRangeUser( minY_, maxY_ );

    hs->Draw("HISTSAME");
    histInfo_data_.h->Draw(drawOp);
    h_format->Draw("AXISSAME");
    legend->Draw();

    DrawLatexAll();

    // -- bottom pad
    c_->cd();
    bottomPad_->cd();

    // -- setup ratio (data/MC) histogram
    SetRatioHistogram();

    Int_t colorRatio = histInfo_data_.color; // -- same color with the data

    h_ratio_dataToStack_->Draw(drawOp);
    h_ratio_dataToStack_->SetStats(kFALSE);
    h_ratio_dataToStack_->SetMarkerStyle(20);
    h_ratio_dataToStack_->SetMarkerColor(colorRatio);
    h_ratio_dataToStack_->SetLineColor(colorRatio);
    h_ratio_dataToStack_->SetFillColorAlpha(kWhite, 0); 
    h_ratio_dataToStack_->SetTitle("");
    PlotTool::SetAxis_BottomPad(h_ratio_dataToStack_->GetXaxis(), h_ratio_dataToStack_->GetYaxis(), titleX_, titleRatio_);
    if( setRangeRatio_ ) h_ratio_dataToStack_->GetYaxis()->SetRangeUser( minRatio_, maxRatio_ );

    TF1 *f_line;
    PlotTool::DrawLine(f_line);

    c_->SaveAs(".pdf");
  }


private:
  void SetMCStack(TLegend *legend)
  {
    hs = new THStack("hs", "");

    Int_t nHistStack = (Int_t)histInfos_.size();
    for(Int_t i=0; i<nHistStack; i++)
    {
      TH1D*& h    = histInfos_[i].h;
      Int_t color = histInfos_[i].color;

      if( setRebin_ ) h->Rebin( nRebin_ );

      h->SetStats(kFALSE);
      h->SetMarkerStyle(20);
      h->SetMarkerColor(color);
      h->SetLineColor(color);
      h->SetFillColor(color);
      h->SetTitle("");

      if( setRangeX_ ) h->GetXaxis()->SetRangeUser( minX_, maxX_ );
      if( setRangeY_ ) h->GetYaxis()->SetRangeUser( minY_, maxY_ );

      hs->Add( h );
    }

    for(Int_t i=nHistStack-1; i>=0; i--) // -- reverse order
      legend->AddEntry(histInfos_[i].h, histInfos_[i].legend);
  }

  void SetDataHistogram(TLegend *legend)
  {
    TH1D*& h           = histInfo_data_.h;
    TString legendName = histInfo_data_.legend;
    Int_t color        = histInfo_data_.color;

    if( setRebin_ ) h->Rebin( nRebin_ );

    h->SetStats(kFALSE);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColorAlpha(kWhite, 0); 
    h->SetTitle("");

    if( setRangeX_ ) h->GetXaxis()->SetRangeUser( minX_, maxX_ );
    if( setRangeY_ ) h->GetYaxis()->SetRangeUser( minY_, maxY_ );

    legend->AddEntry( h, legendName, "EP" ); // -- no horizontal error bar
  }

  void SetRatioHistogram()
  {
    TH1D* h_data = (TH1D*)histInfo_data_.h->Clone();
    h_data->Sumw2();

    TH1D *h_totStack = NULL;
    Int_t nHistStack = (Int_t)histInfos_.size();
    for(Int_t i_stack=0; i_stack<nHistStack; i_stack++)
    {
      histInfos_[i_stack].h->Sumw2();

      if( h_totStack == NULL )
        h_totStack = (TH1D*)histInfos_[i_stack].h->Clone();
      else
        h_totStack->Add( histInfos_[i_stack].h );
    }

    h_ratio_dataToStack_ = (TH1D*)h_data->Clone();
    h_ratio_dataToStack_->Divide( h_data, h_totStack );
  }

};

class GraphCanvas: public CanvasBase
{
public:
  vector<GraphInfo> graphInfos_;

  GraphCanvas()
  {
    // -- member variables are initialized by Init() in CanvasBase()
  }

  GraphCanvas(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE ): GraphCanvas()
  {
    canvasName_ = canvasName;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void Register( TGraphAsymmErrors* g, TString legend, Int_t color  )
  {
    GraphInfo graphInfo{ (TGraphAsymmErrors*)g->Clone(), legend, color };
    graphInfos_.push_back( graphInfo );
  }

  void Draw( TString drawOp = "EPSAME" )
  {
    if( !drawOp.Contains("SAME") ) drawOp = drawOp + "SAME";

    TLegend *legend;
    PlotTool::SetLegend( legend, legendMinX_, legendMinY_, legendMaxX_, legendMaxY_ );
    if( setLegendColumn_ ) legend->SetNColumns(nLegendColumn_);

    // -- draw canvas
    SetCanvas_Square();

    c_->cd();

    Int_t nGraph = graphInfos_.size();
    for(Int_t i=0; i<nGraph; i++)
    {
      TGraphAsymmErrors*& g = graphInfos_[i].g;
      TString legendName = graphInfos_[i].legend;
      Int_t color = graphInfos_[i].color;

      if( i == 0) g->Draw("A"+drawOp);
      else        g->Draw(drawOp);

      g->SetMarkerStyle(20);
      g->SetMarkerColor(color);
      g->SetMarkerSize(1.3);

      g->SetLineColor(color);
      g->SetLineWidth(1.0);

      g->SetFillColorAlpha(kWhite, 0); 
      g->SetTitle("");

      if( i == 0 ) PlotTool::SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX_, titleY_ );
      if( setRangeX_ ) g->GetXaxis()->SetLimits( minX_, maxX_ );
      if( setRangeY_ ) g->GetYaxis()->SetRangeUser( minY_, maxY_ );

      legend->AddEntry( g, legendName );
    }

    legend->Draw();

    DrawLatexAll();

    c_->SaveAs(".pdf");
  }
};

class GraphCanvaswRatio: public GraphCanvas
{
public:
  vector<GraphInfo> graphInfoRatios_;

  GraphCanvaswRatio()
  {
    // -- member variables are initialized by Init() in HistCanvasBase()
  }

  GraphCanvaswRatio(TString canvasName, Bool_t isLogX = kFALSE, Bool_t isLogY = kFALSE ): GraphCanvaswRatio()
  {
    canvasName_ = canvasName;
    isLogX_ = isLogX;
    isLogY_ = isLogY;
  }

  void Draw( TString drawOp = "EPSAME" )
  {
    if( !drawOp.Contains("SAME") ) drawOp = drawOp + "SAME";

    TLegend *legend;
    PlotTool::SetLegend( legend, legendMinX_, legendMinY_, legendMaxX_, legendMaxY_ );
    if( setLegendColumn_ ) legend->SetNColumns(nLegendColumn_);

    // -- draw canvas
    SetCanvas_Ratio();

    c_->cd();
    topPad_->cd();

    Int_t nGraph = graphInfos_.size();
    for(Int_t i=0; i<nGraph; i++)
    {
      TGraphAsymmErrors*& g = graphInfos_[i].g;
      TString legendName = graphInfos_[i].legend;
      Int_t color = graphInfos_[i].color;

      if( i == 0) g->Draw("A"+drawOp);
      else        g->Draw(drawOp);

      g->SetMarkerStyle(20);
      g->SetMarkerColor(color);
      g->SetMarkerSize(1.3);

      g->SetLineColor(color);
      g->SetLineWidth(1.0);

      g->SetFillColorAlpha(kWhite, 0); 
      g->SetTitle("");

      if( i == 0 ) PlotTool::SetAxis_TopPad( g->GetXaxis(), g->GetYaxis(), titleY_ );
      if( setRangeX_ ) g->GetXaxis()->SetLimits( minX_, maxX_ );
      if( setRangeY_ ) g->GetYaxis()->SetRangeUser( minY_, maxY_ );

      legend->AddEntry( g, legendName );
    }

    legend->Draw();

    DrawLatexAll();

    // -- bottom pad
    c_->cd();
    bottomPad_->cd();

    CalcRatioGraph();

    Int_t nGraphRatio = graphInfoRatios_.size();
    for(Int_t i=0; i<nGraphRatio; i++)
    {
      TGraphAsymmErrors*& g_ratio = graphInfoRatios_[i].g;
      Int_t               color   = graphInfoRatios_[i].color;

      if( i == 0) g_ratio->Draw("A"+drawOp);
      else        g_ratio->Draw(drawOp);

      g_ratio->SetMarkerStyle(20);
      g_ratio->SetMarkerColor(color);
      g_ratio->SetMarkerSize(1.3);

      g_ratio->SetLineColor(color);
      g_ratio->SetLineWidth(1.0);

      g_ratio->SetFillColorAlpha(kWhite, 0); 
      g_ratio->SetTitle("");

      if( i == 0 ) SetAxis_BottomPad(g_ratio->GetXaxis(), g_ratio->GetYaxis(), titleX_, titleRatio_);
      if( setRangeX_ )     g_ratio->GetXaxis()->SetLimits( minX_, maxX_ );
      if( setRangeRatio_ ) g_ratio->GetYaxis()->SetRangeUser( minRatio_, maxRatio_ );
    }

    // -- remove points after drawing all of them
    if( doRemoveZeroPoint_ )
    {
      for(Int_t i=0; i<nGraph; i++)
      {
        TGraphAsymmErrors*& g = graphInfos_[i].g;
        RemoveZeroPointInGraph(g);
      }
      for(Int_t i=0; i<nGraphRatio; i++)
      {
        TGraphAsymmErrors*& g_ratio = graphInfoRatios_[i].g;
        RemoveZeroPointInGraph(g_ratio);
      }
    }

    TF1 *f_line;
    PlotTool::DrawLine(f_line);

    c_->SaveAs(".pdf");
  }

  void CalcRatioGraph()
  {
    TGraphAsymmErrors* g_ref = graphInfos_[0].g;

    Int_t nGraph = graphInfos_.size();
    for(Int_t i=1; i<nGraph; i++) // -- starts with 1 -- //
    {
      TGraphAsymmErrors* g_target = (TGraphAsymmErrors*)graphInfos_[i].g->Clone();

      TString legend = graphInfos_[i].legend;
      Int_t color = graphInfos_[i].color;

      TGraphAsymmErrors *g_ratioTemp = MakeRatioGraph( g_target, g_ref );
      if( doRemoveZeroPoint_ ) RemoveZeroPointInGraph(g_ratioTemp);

      GraphInfo graphInfoRatio{ g_ratioTemp, legend, color };
      graphInfoRatios_.push_back( graphInfoRatio );
    }
  }

  // -- NUM = Numerator
  // -- DEN = Denominator
  TGraphAsymmErrors* MakeRatioGraph(TGraphAsymmErrors *g_NUM, TGraphAsymmErrors *g_DEN)
  {
    TGraphAsymmErrors* g_ratio = (TGraphAsymmErrors*)g_DEN->Clone();
    g_ratio->Set(0); // --Remove all points (reset) -- //

    Int_t nPoint_NUM = g_NUM->GetN();
    Int_t nPoint_DEN = g_DEN->GetN();
    if( nPoint_NUM != nPoint_DEN )
      printf("# points is different bewteen two graph...be careful for the ratio plot\n");

    for(Int_t i_p=0; i_p<nPoint_NUM; i_p++)
    {
      Double_t x_NUM, y_NUM;
      g_NUM->GetPoint(i_p, x_NUM, y_NUM);
      Double_t error_NUMLow = g_NUM->GetErrorYlow(i_p);
      Double_t error_NUMHigh = g_NUM->GetErrorYhigh(i_p);
      // -- take the larger uncertainty
      Double_t error_NUM = error_NUMLow > error_NUMHigh ? error_NUMLow : error_NUMHigh;

      Double_t x_DEN, y_DEN;
      g_DEN->GetPoint(i_p, x_DEN, y_DEN);
      Double_t error_DENLow = g_DEN->GetErrorYlow(i_p);
      Double_t error_DENHigh = g_DEN->GetErrorYhigh(i_p);
      // -- take the larger uncertainty
      Double_t error_DEN = error_DENLow > error_DENHigh ? error_DENLow : error_DENHigh;

      Double_t ratio;
      Double_t ratio_error;
      if( (nPoint_NUM != nPoint_DEN) && i_p >= nPoint_DEN )
      {
        ratio = 0;
        ratio_error = 0;
      }
      // else if(y_Type1 != 0 && error_Type1 != 0 && y_Type2 != 0 && error_Type2 != 0)
      else if(y_DEN != 0)
      {
        ratio = y_NUM / y_DEN;
        ratio_error = this->Error_PropagatedAoverB(y_NUM, error_NUM, y_DEN, error_DEN);
        //calculate Scale Factor(Type1/Type2) & error

        // cout << "ratio: " << ratio << " ratio_error: " << ratio_error << endl;
      }
      else
      {
        cout << "Denominator is 0! ... ratio and its error are set as 0" << endl;
        ratio = 0;
        ratio_error = 0;
      }

      //Set Central value
      g_ratio->SetPoint(i_p, x_NUM, ratio);

      //Set the error
      Double_t error_XLow = g_NUM->GetErrorXlow(i_p);
      Double_t error_Xhigh = g_NUM->GetErrorXhigh(i_p);
      g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);

      // cout << endl;
    }

    return g_ratio;
  }

  Double_t Error_PropagatedAoverB(Double_t A, Double_t sigma_A, Double_t B, Double_t sigma_B)
  {
    Double_t ratio_A = (sigma_A) / A;
    Double_t ratio_B = (sigma_B) / B;

    Double_t errorSquare = ratio_A * ratio_A + ratio_B * ratio_B;

    return (A/B) * sqrt(errorSquare);
  }
};

}; // -- namespace PlotTool