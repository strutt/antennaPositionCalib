{

  gStyle->SetOptFit(01111);

  auto cv = new TChain("angResTree");
  auto ch = new TChain("angResTree");

  // cv->Add("vpolLindasHpolValues/*");
  cv->Add("vpolLindasHpolValuesFull/*");  
  ch->Add("hpolLindasHpolValues/*");


  auto c1 = new TCanvas();  
  c1->SetLogy(1);

  ch->Draw("deltaThetaDeg>>hThetaHPOL", "TMath::Abs(deltaThetaDeg) < 5", "");  
  auto hThetaHPOL = (TH1F*) gDirectory->Get("hThetaHPOL");
  cv->Draw("deltaThetaDeg>>hThetaVPOL", "TMath::Abs(deltaThetaDeg) < 5", "sames");
  auto hThetaVPOL = (TH1F*) gDirectory->Get("hThetaVPOL");

  auto c2 = new TCanvas();
  c2->SetLogy(1);

  ch->Draw("deltaPhiDeg>>hPhiHPOL", "TMath::Abs(deltaPhiDeg) < 5", "");  
  auto hPhiHPOL = (TH1F*) gDirectory->Get("hPhiHPOL");
  cv->Draw("deltaPhiDeg>>hPhiVPOL", "TMath::Abs(deltaPhiDeg) < 5", "sames");
  auto hPhiVPOL = (TH1F*) gDirectory->Get("hPhiVPOL");


  std::cout << hThetaHPOL->Integral() << std::endl;
  std::cout << hPhiHPOL->Integral() << std::endl;    
  std::cout << hThetaVPOL->Integral() << std::endl;
  std::cout << hPhiVPOL->Integral() << std::endl;  

  
  hThetaHPOL->SetLineColor(kRed);
  
  c1->cd();
  
  hThetaHPOL->SetLineColor(kRed);
  hThetaHPOL->Scale(1./hThetaHPOL->Integral());
  hThetaHPOL->SetTitle("#delta#theta for VPOL (LDB) and HPOL (WAIS) pulses; #theta_{Expected} - #theta_{Measured} (Degrees); Events per bin (normalized)");
  hThetaHPOL->Fit("gaus");
  auto fThetaHPOL = (TF1*) hThetaHPOL->FindObject("gaus");
  fThetaHPOL->SetName("fThetaHPOL");
  fThetaHPOL->SetLineColor(kRed);
  
  
  hThetaVPOL->SetLineColor(kBlue);  
  hThetaVPOL->Scale(1./hThetaVPOL->Integral());  
  hThetaVPOL->Fit("gaus");
  auto fThetaVPOL = (TF1*) hThetaVPOL->FindObject("gaus");
  fThetaVPOL->SetName("fThetaVPOL");
  fThetaVPOL->SetLineColor(kBlue);
  c1->Update();
  
  auto lTheta = new TLegend(0.8, 0.8, 1, 1);
  lTheta->AddEntry(hThetaHPOL, "#delta#theta HPOL", "l");
  lTheta->AddEntry(hThetaVPOL, "#delta#theta VPOL", "l");  
  lTheta->Draw();
  
  auto sThetaVPOL = (TPaveStats*) hThetaVPOL->FindObject("stats");
  sThetaVPOL->SetTextColor(kBlue);
  sThetaVPOL->SetX1NDC(0.8);
  sThetaVPOL->SetX2NDC(1.0);    
  sThetaVPOL->SetY1NDC(0.2);
  sThetaVPOL->SetY2NDC(0.5);

  auto sThetaHPOL = (TPaveStats*) hThetaHPOL->FindObject("stats");
  sThetaHPOL->SetTextColor(kRed);
  sThetaHPOL->SetX1NDC(0.8);
  sThetaHPOL->SetX2NDC(1.0);    
  sThetaHPOL->SetY1NDC(0.5);
  sThetaHPOL->SetY2NDC(0.8);




  c2->cd();
  
  hPhiHPOL->SetLineColor(kRed);
  hPhiHPOL->Scale(1./hPhiHPOL->Integral());
  hPhiHPOL->SetTitle("#delta#phi for VPOL (LDB) and HPOL (WAIS) pulses; #phi_{Expected} - #phi_{Measured} (Degrees); Events per bin (normalized)");
  hPhiHPOL->Fit("gaus");
  auto fPhiHPOL = (TF1*) hPhiHPOL->FindObject("gaus");
  fPhiHPOL->SetName("fPhiHPOL");
  fPhiHPOL->SetLineColor(kRed);
  
  
  hPhiVPOL->SetLineColor(kBlue);  
  hPhiVPOL->Scale(1./hPhiVPOL->Integral());  
  hPhiVPOL->Fit("gaus");
  auto fPhiVPOL = (TF1*) hPhiVPOL->FindObject("gaus");
  fPhiVPOL->SetName("fPhiVPOL");
  fPhiVPOL->SetLineColor(kBlue);
  c2->Update();

  
  auto lPhi = new TLegend(0.8, 0.8, 1, 1);
  lPhi->AddEntry(hPhiHPOL, "#delta#phi HPOL", "l");
  lPhi->AddEntry(hPhiVPOL, "#delta#phi VPOL", "l");  
  lPhi->Draw();
  
  auto sPhiVPOL = (TPaveStats*) hPhiVPOL->FindObject("stats");
  sPhiVPOL->SetTextColor(kBlue);
  sPhiVPOL->SetX1NDC(0.8);
  sPhiVPOL->SetX2NDC(1.0);    
  sPhiVPOL->SetY1NDC(0.2);
  sPhiVPOL->SetY2NDC(0.5);

  auto sPhiHPOL = (TPaveStats*) hPhiHPOL->FindObject("stats");
  sPhiHPOL->SetTextColor(kRed);
  sPhiHPOL->SetX1NDC(0.8);
  sPhiHPOL->SetX2NDC(1.0);    
  sPhiHPOL->SetY1NDC(0.5);
  sPhiHPOL->SetY2NDC(0.8);



  
  auto c3 = new TCanvas();
  c3->Divide(2, 2);
  c3->cd(1);
  gPad->SetLogz(1);
  ch->Draw("deltaThetaDeg:thetaExpected>>hThetaHPOL2", "TMath::Abs(deltaThetaDeg) < 3", "colz");  
  auto hThetaHPOL2 = (TH1F*) gDirectory->Get("hThetaHPOL2");
  hThetaHPOL2->SetTitle("#delta#theta for HPOL (WAIS) pulses; #theta_{Expected} (Degrees); #theta_{Expected} - #theta_{Measured} (Degrees)");
  c3->cd(2);
  gPad->SetLogz(1);
  cv->Draw("deltaThetaDeg:thetaExpected>>hThetaVPOL2", "TMath::Abs(deltaThetaDeg) < 3", "colz");
  auto hThetaVPOL2 = (TH1F*) gDirectory->Get("hThetaVPOL2");
  hThetaVPOL2->SetTitle("#delta#theta for VPOL (LDB) pulses; #theta_{Expected} (Degrees); #theta_{Expected} - #theta_{Measured} (Degrees)");
  
  c3->cd(3);
  gPad->SetLogz(1);
  ch->Draw("deltaPhiDeg:phiExpected>>hPhiHPOL2", "TMath::Abs(deltaPhiDeg) < 3", "colz");  
  auto hPhiHPOL2 = (TH1F*) gDirectory->Get("hPhiHPOL2");
  hPhiHPOL2->SetTitle("#delta#phi for HPOL (WAIS) pulses; #phi_{Expected} (Degrees); #phi_{Expected} - #phi_{Measured} (Degrees)");

  c3->cd(4);
  gPad->SetLogz(1);
  cv->Draw("deltaPhiDeg:phiExpected>>hPhiVPOL2", "TMath::Abs(deltaPhiDeg) < 3", "colz");
  auto hPhiVPOL2 = (TH1F*) gDirectory->Get("hPhiVPOL2");
  hPhiVPOL2->SetTitle("#delta#phi for VPOL (LDB) pulses; #phi_{Expected} (Degrees); #phi_{Expected} - #phi_{Measured} (Degrees)");  

  
}
