void drawAngularResolutionPlots2(){

  gStyle->SetOptFit(01111);
  
  Int_t colCorrected = 2;
  Int_t colNotCorrected = 4;  
  Bool_t friends = kTRUE;
  auto chain = new TChain("angResTree");
  auto chain2 = new TChain("angResTreeWithDynamicPitchRoll");  

  // chain->Add("generateAngularResolutionTree_run*.root");
  // chain->Add("generateAngularResolutionTree_run342*.root");
  // chain->Add("generateAngularResolutionTree_run352*Plots.root");
  // chain->Add("generateAngularResolutionTree_run*Plots_noPitchRoll.root");
  // chain->Add("generateAngularResolutionTree_run342-342Plots_part*.root");
  for(int run=331; run<=354; run++){
    TString fileName = TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    chain->Add(fileName);
    fileName = TString::Format("generateAngularResolutionTreeFriend_run%d-%dPlots.root", run, run);    
    chain2->Add(fileName);
  }
  chain->AddFriend(chain2);
  
  //WithProvisionalBestConstantOffsets
  //WithProvisionalBestConstantOffsets
  
  auto c1 = new TCanvas();

  chain->Draw("deltaPhiDegWithDynamicPitchRoll", "TMath::Abs(deltaPhiDegWithDynamicPitchRoll) < 3", "goff");
  TH1F* hDeltaPhiDeg3 = (TH1F*) gROOT->Get("htemp")->Clone("hDeltaPhiDegWithDynamicPitchRoll");
  hDeltaPhiDeg3->SetTitle("Close up of #phi_{expected} - #phi_{zoom}; #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg3->Draw(); 
  hDeltaPhiDeg3->GetYaxis()->SetNoExponent(1);
  hDeltaPhiDeg3->Fit("gaus", "0");
  
  chain->Draw("deltaPhiDegWithProvisionalBestConstantOffsets", "TMath::Abs(deltaPhiDegWithProvisionalBestConstantOffsets) < 3", "goff");
  TH1F* hDeltaPhiDeg = (TH1F*) gROOT->Get("htemp")->Clone("hDeltaPhiDeg");
  hDeltaPhiDeg->SetTitle("Close up of #phi_{expected} - #phi_{zoom}; #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg->Draw("sames"); 
  hDeltaPhiDeg->GetYaxis()->SetNoExponent(1);
  hDeltaPhiDeg->Fit("gaus");

  
  hDeltaPhiDeg->SetLineColor(colNotCorrected);
  hDeltaPhiDeg->Draw();
  hDeltaPhiDeg->GetYaxis()->SetNoExponent(1);
  hDeltaPhiDeg->Fit("gaus", "0");
  hDeltaPhiDeg->GetFunction("gaus")->SetLineColor(hDeltaPhiDeg->GetLineColor());
  hDeltaPhiDeg->GetFunction("gaus")->SetLineStyle(2);
  hDeltaPhiDeg->GetFunction("gaus")->Draw("sames");
  hDeltaPhiDeg3->SetLineColor(colCorrected);
  hDeltaPhiDeg3->Draw("sames");
  hDeltaPhiDeg3->GetYaxis()->SetNoExponent(1);
  hDeltaPhiDeg3->Fit("gaus", "0");
  hDeltaPhiDeg3->GetFunction("gaus")->SetLineColor(hDeltaPhiDeg3->GetLineColor());
  hDeltaPhiDeg3->GetFunction("gaus")->SetLineStyle(2);
  hDeltaPhiDeg3->GetFunction("gaus")->Draw("sames");


  c1->Update();

  auto s1a = (TPaveStats*) hDeltaPhiDeg3->FindObject("stats");
  s1a->SetTextColor(hDeltaPhiDeg3->GetLineColor());
  s1a->SetX1NDC(0.6);
  s1a->SetY1NDC(0.6);
  s1a->SetX2NDC(0.98);
  s1a->SetY2NDC(0.8);  

  auto s2a = (TPaveStats*) hDeltaPhiDeg->FindObject("stats");
  s2a->SetTextColor(hDeltaPhiDeg->GetLineColor());
  s2a->SetX1NDC(0.6);
  s2a->SetY1NDC(0.4);
  s2a->SetX2NDC(0.98);
  s2a->SetY2NDC(0.6);  

  auto l2 = new TLegend(0.6, 0.8, 0.98, 1);
  l2->AddEntry(hDeltaPhiDeg3, "With dynamic pitch/roll", "l");
  l2->AddEntry(hDeltaPhiDeg, "With best constant offsets", "l");
  l2->SetLineWidth(0);  
  l2->Draw();
  
  hDeltaPhiDeg->GetXaxis()->SetRangeUser(-4, 4);  
  c1->Modified();
  c1->Update();
    

  

  
  // auto c0 = new TCanvas();
  // chain->Draw("deltaPhiDegWithDynamicPitchRoll:zoomPhiDeg", "TMath::Abs(deltaPhiDeg) < 3", "colz");

  // return;

  auto c2 = new TCanvas();  
  chain->Draw("deltaPhiDegWithProvisionalBestConstantOffsets:zoomPhiDeg", "TMath::Abs(deltaPhiDegWithProvisionalBestConstantOffsets) < 3", "goff");  
  TH2F* hDeltaPhiDeg_vs_phiZoom = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDeg_vs_phiZoom");
  hDeltaPhiDeg_vs_phiZoom->SetTitle("Close up of #phi_{expected} - #phi_{zoom}; #phi_{zoom} (Degrees); #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg_vs_phiZoom->Draw("colz");

  TProfile* hDeltaPhiDeg_vs_phiZoom_pfx = hDeltaPhiDeg_vs_phiZoom->ProfileX();
  hDeltaPhiDeg_vs_phiZoom_pfx->Draw();
  hDeltaPhiDeg_vs_phiZoom_pfx->SetLineColor(colNotCorrected);
  hDeltaPhiDeg_vs_phiZoom_pfx->GetYaxis()->SetRangeUser(-2, 2);  
  
  chain->Draw("deltaPhiDegWithDynamicPitchRoll:zoomPhiDeg", "TMath::Abs(deltaPhiDegWithDynamicPitchRoll) < 3", "goff");
  TH2F* hDeltaPhiDeg_vs_phiZoom3 = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDegWithDynamicPitchRoll_vs_phiZoom");
  TProfile* hDeltaPhiDeg_vs_phiZoom3_pfx = hDeltaPhiDeg_vs_phiZoom3->ProfileX();
  hDeltaPhiDeg_vs_phiZoom3_pfx->Draw("same");
  hDeltaPhiDeg_vs_phiZoom3_pfx->SetLineColor(colCorrected);
  l2->Draw();
  
  

  auto c3 = new TCanvas();
  chain->Draw("deltaThetaDegWithProvisionalBestConstantOffsets", "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 3", "goff");
  TH1F* hDeltaThetaDeg = (TH1F*) gROOT->Get("htemp")->Clone("hDeltaThetaDeg");
  hDeltaThetaDeg->SetTitle("Close up of #theta_{expected} - #theta_{zoom}; #delta#theta (degrees); Events per bin");
  hDeltaThetaDeg->Rebin(2);

  chain->Draw("deltaThetaDegWithDynamicPitchRoll", "TMath::Abs(deltaThetaDegWithDynamicPitchRoll) < 3", "goff");
  TH1F* hDeltaThetaDeg3 = (TH1F*) gROOT->Get("htemp")->Clone("hDeltaThetaDegWithDynamicPitchRoll");
  hDeltaThetaDeg3->SetTitle("Close up of #theta_{expected} - #theta_{zoom} Best fit offsets; #delta#theta (degrees); Events per bin");
  hDeltaThetaDeg3->Rebin(2);
  
  hDeltaThetaDeg->SetLineColor(colNotCorrected);
  hDeltaThetaDeg->Draw();
  hDeltaThetaDeg->GetYaxis()->SetNoExponent(1);
  hDeltaThetaDeg->Fit("gaus", "0");
  hDeltaThetaDeg->GetFunction("gaus")->SetLineColor(hDeltaThetaDeg->GetLineColor());  
  hDeltaThetaDeg->GetFunction("gaus")->SetLineStyle(2);
  hDeltaThetaDeg->GetFunction("gaus")->Draw("sames");  
  hDeltaThetaDeg3->SetLineColor(colCorrected);
  hDeltaThetaDeg3->Draw("sames");
  hDeltaThetaDeg3->GetYaxis()->SetNoExponent(1);
  hDeltaThetaDeg3->Fit("gaus", "0");
  hDeltaThetaDeg3->GetFunction("gaus")->SetLineColor(hDeltaThetaDeg3->GetLineColor());
  hDeltaThetaDeg3->GetFunction("gaus")->SetLineStyle(2);
  hDeltaThetaDeg3->GetFunction("gaus")->Draw("sames");
  c3->Update();
  

  auto s1 = (TPaveStats*) hDeltaThetaDeg3->FindObject("stats");
  s1->SetTextColor(hDeltaThetaDeg3->GetLineColor());
  s1->SetX1NDC(0.6);
  s1->SetY1NDC(0.6);
  s1->SetX2NDC(0.98);
  s1->SetY2NDC(0.8);  

  auto s2 = (TPaveStats*) hDeltaThetaDeg->FindObject("stats");
  s2->SetTextColor(hDeltaThetaDeg->GetLineColor());
  s2->SetX1NDC(0.6);
  s2->SetY1NDC(0.4);
  s2->SetX2NDC(0.98);
  s2->SetY2NDC(0.6);  

  auto l3 = new TLegend(0.6, 0.8, 0.98, 1);
  l3->AddEntry(hDeltaThetaDeg3, "With dynamic pitch/roll", "l");
  l3->AddEntry(hDeltaThetaDeg, "With best constant offsets", "l");
  l3->SetLineWidth(0);  
  l3->Draw();
  
  hDeltaThetaDeg->GetXaxis()->SetRangeUser(-2, 2);  
  c3->Modified();
  c3->Update();
  
  
  auto c4 = new TCanvas();
  // chain->Draw("deltaThetaDegWithProvisionalBestConstantOffsets:zoomPhiDeg", "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 3", "goff");
  chain->Draw("deltaThetaDegWithProvisionalBestConstantOffsets:zoomPhiDeg", "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 3", "goff"); 
  TH2F* hDeltaThetaDeg_vs_phiZoom = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDeg_vs_phiZoom");
  hDeltaThetaDeg_vs_phiZoom->SetTitle("Close up of #theta_{expected} - #theta_{zoom}; #phi_{zoom} (Degrees); #delta#theta (degrees); Events per bin");
  TProfile* hDeltaThetaDeg_vs_phiZoom_pfx = hDeltaThetaDeg_vs_phiZoom->ProfileX();

  chain->Draw("deltaThetaDegWithDynamicPitchRoll:zoomPhiDeg", "TMath::Abs(deltaThetaDegWithDynamicPitchRoll) < 3", "goff"); 
  TH2F* hDeltaThetaDeg_vs_phiZoom3 = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDegWithDynamicPitchRoll_vs_phiZoom");
  hDeltaThetaDeg_vs_phiZoom3->SetTitle("Close up of #theta_{expected} - #theta_{zoom}; #phi_{zoom} (Degrees); #delta#theta (degrees); Events per bin");
  hDeltaThetaDeg_vs_phiZoom3->Draw("colz");

  TProfile* hDeltaThetaDeg_vs_phiZoom3_pfx = hDeltaThetaDeg_vs_phiZoom3->ProfileX();
  hDeltaThetaDeg_vs_phiZoom_pfx->SetLineColor(colNotCorrected);
  hDeltaThetaDeg_vs_phiZoom_pfx->Draw();
  hDeltaThetaDeg_vs_phiZoom3_pfx->SetLineColor(colCorrected);
  hDeltaThetaDeg_vs_phiZoom3_pfx->Draw("same");
  hDeltaThetaDeg_vs_phiZoom_pfx->GetYaxis()->SetRangeUser(-0.5, 0.5);
  l2->Draw();
}
