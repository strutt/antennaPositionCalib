{

  auto cBefore = new TChain("angResTree");
  auto cAfter = new TChain("angResTree");

  cBefore->Add("photogrammetryNoPitchRoll/generateAngularResolutionTreePlots_*_2015-10-26_11-23*");
  cAfter->Add("hpolLindasHpolValues/generateAngularResolutionTreePlots_*_2015-10-26_13*");  
  
  // cAfter->Draw("deltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
  cAfter->Draw("zoomPeak", "TMath::Abs(deltaThetaDeg) < 3", "goff");  
  auto hThetaAfter = (TH1F*) gROOT->Get("htemp")->Clone("hThetaAfter");
  hThetaAfter->SetLineColor(kBlue);
  hThetaAfter->GetYaxis()->SetTitle("Entries per bin");
  hThetaAfter->GetXaxis()->SetTitle("#delta#theta (degrees)");  
    
  // cBefore->Draw("deltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
  cBefore->Draw("zoomPeak", "TMath::Abs(deltaThetaDeg) < 3", "goff");  
  auto hThetaBefore = (TH1F*) gROOT->Get("htemp")->Clone("hThetaBefore");
  hThetaBefore->SetLineColor(kRed);


  auto c1 = new TCanvas();

  hThetaAfter->Draw();
  hThetaBefore->Draw("same");    

  auto l = new TLegend(0.7, 0.8, 0.99, 0.99);
  l->AddEntry(hThetaAfter, "Linda Phase Centers + extra delays", "l");
  l->AddEntry(hThetaBefore, "Photogrammetry", "l");  
  l->Draw();

  cAfter->Draw("deltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 3", "goff");
  auto hPhiAfter = (TH1F*) gROOT->Get("htemp")->Clone("hPhiAfter");
  hPhiAfter->SetLineColor(kBlue);
  hPhiAfter->GetYaxis()->SetTitle("Entries per bin");
  hPhiAfter->GetXaxis()->SetTitle("#delta#phi (degrees)");  

  
  cBefore->Draw("deltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 3", "goff");
  auto hPhiBefore = (TH1F*) gROOT->Get("htemp")->Clone("hPhiBefore");
  hPhiBefore->SetLineColor(kRed);

  
  auto c2 = new TCanvas();
  
  hPhiAfter->Draw();
  hPhiBefore->Draw("same");  
  l->Draw();






}
