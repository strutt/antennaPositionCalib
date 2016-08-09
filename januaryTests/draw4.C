void draw4(){

  gStyle->SetOptStat("mre");


  // TFile* f = TFile::Open("testing/corInterp/latestNums/generateAngularResolutionTreePlots_352_2016-02-21_15-38-05.root"); //generateAngularResolutionTreePlots_352_2016-02-19_18-01-52.root");
  TFile* f = TFile::Open("../generateAngularResolutionTreePlots_352_2016-02-22_14-07-28.root");
  
  auto angResTree = (TTree*) f->Get("angResTree");

  const int numHists = 3;
  std::vector<TH1D*> hists;

  const Double_t maxAngle = 5;
  const Int_t numBins = 64;
  TH1D* h1 = new TH1D("h1", "globalPhiDeg", numBins, -maxAngle, maxAngle);
  TH1D* h2 = new TH1D("h2", "triggeredPhiDeg", numBins, -maxAngle, maxAngle);
  TH1D* h3 = new TH1D("h3", "zoomPhiDeg", numBins, -maxAngle, maxAngle);

  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h3->SetLineColor(kMagenta);  
  
  hists.push_back(h1);
  hists.push_back(h2);
  hists.push_back(h3);

  angResTree->Draw("phiExpected - globalPhiDeg>>h1", "", "goff");
  angResTree->Draw("phiExpected - triggeredPhiDeg>>h2", "", "goff");
  angResTree->Draw("deltaPhiDeg>>h3", "", "goff");

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());  
  
  TCanvas* c1 = RootTools::drawHistsWithStatsBoxes(numHists, &hists[0],
						   "", "mre");
  auto l = c1->BuildLegend();
  // h1->SetMaximum(0.25);
  h1->SetMaximum(0.12);  

  h1->SetTitle("Comparison of different my different reconstruction options (run 352); #delta#phi (Degrees); Fraction of events;");
  
}
