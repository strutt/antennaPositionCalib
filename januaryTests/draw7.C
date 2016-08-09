void draw7(){

  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_11-57-40.root");
  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_13-05-06.root");
  // TTree* t = (TTree*) f->Get("angResTree");

  TChain* t = new TChain("angResTree");
  t->Add("testing/comboStudies/generateAngularResolutionTreePlots_*");

  
  const int numBins = 64;
  TH2D* h2 = new TH2D("h2", "Option 3", numBins, 6, 7.5, numBins, -5, 5);
  // t->Draw("deltaPhiDeg3:zoomThetaDeg3>>h2", "TMath::Abs(deltaPhiDeg3) < 5 && zoomThetaDeg3 < 0", "colz");
  // t->Draw("deltaPhiDeg:zoomThetaDeg>>h2", "TMath::Abs(deltaPhiDeg) < 5 && zoomThetaDeg < 0", "colz");
  // t->Draw("deltaPhiDeg:-1*zoomThetaDeg>>h2", "TMath::Abs(deltaPhiDeg) < 5", "colz");
  t->Draw("deltaPhiDeg:-1*thetaExpected>>h2", "TMath::Abs(deltaPhiDeg) < 5", "colz");

  auto h2_pfx = h2->ProfileX();
  new TCanvas();
  h2_pfx->Draw();

  h2_pfx->SetTitle("#delta#phi (Degrees) as a function of #theta (all WAIS pulses); #theta_{expected} (Degrees) (same way around as Cosmin); #delta#phi (Degrees)");
  // h2_pfx->SetTitle("#delta#phi (Degrees) as a function of #theta (all WAIS pulses); #theta (Degrees); #delta#phi (Degrees)");

  new TCanvas();
  t->Draw("deltaThetaDeg3", "TMath::Abs(deltaThetaDeg3) < 2", "colz");
  // t->Draw("deltaThetaDeg3", "TMath::Abs(deltaThetaDeg3) < 2", "colz");    
  
}
