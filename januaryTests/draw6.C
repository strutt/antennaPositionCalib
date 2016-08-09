void draw6(){

  gStyle->SetOptStat("mre");


  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_11-57-40.root");
  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_13-05-06.root");
  // TTree* t = (TTree*) f->Get("angResTree");

  TChain* t = new TChain("angResTree");
  // t->Add("testing/comboStudies/generateAngularResolutionTreePlots_*");
  t->Add("offAxisDelay/generateAngularResolutionTreePlots_*");  

  

  TH1D* hs[8];
  const int numBins = 128;
  hs[0] = new TH1D("h", "Option 1", numBins, -5, 5);
  t->Draw("deltaPhiDeg>>h", "TMath::Abs(deltaPhiDeg) < 5", "goff");
  hs[0]->SetLineColor(kBlack);
  
  hs[1] = new TH1D("h2", "Option 2", numBins, -5, 5);
  t->Draw("deltaPhiDeg2>>h2", "TMath::Abs(deltaPhiDeg2) < 5", "goff");
  hs[1]->SetLineColor(kRed);
  
  hs[2] = new TH1D("h3", "Option 3", numBins, -5, 5);
  t->Draw("deltaPhiDeg3>>h3", "TMath::Abs(deltaPhiDeg3) < 5", "goff");
  hs[2]->SetLineColor(kBlue);
  
  hs[3] = new TH1D("h4", "Option 4", numBins, -5, 5);
  t->Draw("deltaPhiDeg4>>h4", "TMath::Abs(deltaPhiDeg4) < 5", "goff");
  hs[3]->SetLineColor(kMagenta);
  
  hs[4] = new TH1D("h5", "Option 5", numBins, -5, 5);
  t->Draw("deltaPhiDeg5>>h5", "TMath::Abs(deltaPhiDeg5) < 5", "goff");
  hs[4]->SetLineColor(kOrange+7);
  
  hs[5] = new TH1D("h6", "Option 6", numBins, -5, 5);
  t->Draw("deltaPhiDeg6>>h6", "TMath::Abs(deltaPhiDeg6) < 5", "goff");
  hs[5]->SetLineColor(kGreen+2);

  hs[6] = new TH1D("h7", "Option 7", numBins, -5, 5);
  t->Draw("deltaPhiDeg7>>h7", "TMath::Abs(deltaPhiDeg7) < 5", "goff");
  hs[6]->SetLineColor(kViolet+2);

  hs[7] = new TH1D("h8", "Option 8", numBins, -5, 5);
  t->Draw("deltaPhiDeg8>>h8", "TMath::Abs(deltaPhiDeg8) < 5", "goff");
  hs[7]->SetLineColor(kCyan+2);
    
  TCanvas* c1 = RootTools::drawHistsWithStatsBoxes(8, hs, "", "mre");
  c1->BuildLegend();
  
  // hs[0]->SetTitle("Comparison of different combinatorics (run 352); #delta#phi (Degrees); Events / bin");
  hs[0]->SetTitle("Comparison of different combinatorics (all WAIS runs); #delta#phi (Degrees); Events / bin");  


  auto c2 = new TCanvas();
  // TH2D* h2D = new TH2D("h2D", "Comparison of different combinatorics (run 352); Option; #delta#phi (Degrees); Events / bin", 8, 0, 8, numBins, -5, 5);
  TH2D* h2D = new TH2D("h2D", "Comparison of different combinatorics (all WAIS runs); Option; #delta#phi (Degrees); Events / bin", 8, 0, 8, numBins, -5, 5);  
  for(int optInd=0; optInd<8; optInd++){
    for(int i=1; i<=hs[optInd]->GetNbinsX(); i++){
      h2D->SetBinContent(optInd+1, i, hs[optInd]->GetBinContent(i));
    }
    h2D->GetXaxis()->SetBinLabel(optInd+1, TString::Format("%d", optInd+1));
  }
  h2D->GetXaxis()->SetLabelOffset(0.01);
  h2D->Draw("colz");
}
