void drawThesis(){

  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_11-57-40.root");
  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_13-05-06.root");
  // TTree* t = (TTree*) f->Get("angResTree");

  TChain* t = new TChain("angResTree");
  t->Add("testing/comboStudies/generateAngularResolutionTreePlots_*");

  

  const int numBins = 64; //128;
  const double maxDeg = 2;
  TH1D* hTheta = new TH1D("hTheta", "#delta#theta Distribution of WAIS pulses; #delta#theta (Degrees); Events per bin", numBins, -maxDeg, maxDeg);

  TH1D* hPhi = new TH1D("hPhi", "#delta#phi Distribution of WAIS pulses; #delta#phi (Degrees); Events per bin", numBins, -maxDeg, maxDeg); 

  // t->Draw("deltaPhiDeg>>h2", "TMath::Abs(deltaPhiDeg) < 5", "colz");
  t->Draw("deltaThetaDeg3>>hTheta", "TMath::Abs(deltaPhiDeg3) < 5", "goff");
  t->Draw("deltaPhiDeg3>>hPhi", "TMath::Abs(deltaPhiDeg3) < 5", "goff");

  


  TH1D* hs[2] = {hTheta, hPhi};
  hPhi->SetLineColor(kRed);
  hTheta->SetLineColor(kBlue);  
  auto c1 = RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");

  TLegend* l1 = new TLegend(0.8, 0.8, 1, 1);
  l1->AddEntry(hPhi, "#delta#phi (Degrees)", "l");
  l1->AddEntry(hTheta, "#delta#phi (Degrees)", "l");
  hs[0]->SetTitle("Angular resolution of WAIS pulses");
}
