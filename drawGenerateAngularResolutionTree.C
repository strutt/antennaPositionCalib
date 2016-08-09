{

  auto c = new TChain("angResTree");
  for(int run=331; run<=333; run++){
    auto fName = TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    c->Add(fName);
  }

  TCanvas* c1 = new TCanvas();
  c->Draw("deltaPhiDeg",  "l3TrigPatternH > 0 && heading > -100 && TMath::Abs(deltaPhiDeg) < 2", "colz");
  c1->SetLogy(1);










}
