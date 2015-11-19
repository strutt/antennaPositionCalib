{

  Bool_t friends = kTRUE;
  auto chain = new TChain("angResTree");
  auto chain2 = new TChain("angResTreeWithDynamicPitchRoll");  

  // chain->Add("generateAngularResolutionTree_run*.root");
  // chain->Add("generateAngularResolutionTree_run342*.root");
  // chain->Add("generateAngularResolutionTree_run352*Plots.root");
  // chain->Add("generateAngularResolutionTree_run*Plots_noPitchRoll.root");
  // chain->Add("generateAngularResolutionTree_run342-342Plots_part*.root");
  chain->Add("generateAngularResolutionTree_run*Plots.root");
  chain2->Add("generateAngularResolutionTreeFriend_run*Plots.root");    
  cout << chain->GetEntries() << endl;

  if(friends){
    chain->AddFriend(chain2);
  }

  
  TCanvas* c1 = new TCanvas();

  chain->Draw("deltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 3", "goff");
  TH1F* hDeltaPhiDeg = (TH1F*) htemp->Clone("hDeltaPhiDeg");
  hDeltaPhiDeg->SetTitle("Close up of #phi_{expected} - #phi_{zoom} No dynamic correction; #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg->Draw(); 
  hDeltaPhiDeg->GetYaxis()->SetNoExponent(1);
  hDeltaPhiDeg->Fit("gaus");


  TCanvas* c2 = new TCanvas();
  chain->Draw("deltaPhiDeg:zoomPhiDeg", "TMath::Abs(deltaPhiDeg) < 3", "goff");
  TH2F* hDeltaPhiDeg_vs_phiZoom = (TH2F*) htemp->Clone("hDeltaPhiDeg_vs_phiZoom");
  hDeltaPhiDeg_vs_phiZoom->SetTitle("Close up of #phi_{expected} - #phi_{zoom} No dynamic correction; #phi_{zoom} (Degrees); #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg_vs_phiZoom->Draw("colz");


  TCanvas* c3 = new TCanvas();
  chain->Draw("deltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
  TH1F* hDeltaThetaDeg = (TH1F*) htemp->Clone("hDeltaThetaDeg");
  hDeltaThetaDeg->SetTitle("Close up of #theta_{expected} - #theta_{zoom} No dynamic correction; #delta#theta (degrees); Events per bin");
  hDeltaThetaDeg->Rebin(2);  
  hDeltaThetaDeg->Draw();
  hDeltaThetaDeg->GetYaxis()->SetNoExponent(1);  
  hDeltaThetaDeg->Fit("gaus");  

  
  TCanvas* c4 = new TCanvas();
  // chain->Draw("deltaThetaDeg:zoomPhiDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
  chain->Draw("deltaThetaDeg:zoomPhiDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff"); 
  TH2F* hDeltaThetaDeg_vs_thetaZoom = (TH2F*) htemp->Clone("hDeltaThetaDeg_vs_thetaZoom");
  hDeltaThetaDeg_vs_thetaZoom->SetTitle("Close up of #theta_{expected} - #theta_{zoom} No dynamic correction; #phi_{zoom} (Degrees); #delta#theta (degrees); Events per bin");
  hDeltaThetaDeg_vs_thetaZoom->Draw("colz");



  if(friends){

  //   TCanvas* c1b = new TCanvas();
  //   chain->Draw("deltaPhiDegWithDynamicPitchRoll", "TMath::Abs(deltaPhiDegWithDynamicPitchRoll) < 3", "goff");
  //   TH1F* hDeltaPhiDeg2 = (TH1F*) htemp->Clone("hDeltaPhiDegWithDynamicPitchRoll");
  //   hDeltaPhiDeg2->SetTitle("Close up of #phi_{expected} - #phi_{zoom} With dynamic correction; #delta#phi (degrees); Events per bin");
  //   hDeltaPhiDeg2->Draw(); 
  //   hDeltaPhiDeg2->GetYaxis()->SetNoExponent(1);
  //   hDeltaPhiDeg2->Fit("gaus");

  
  
  // TCanvas* c2b = new TCanvas();
  // chain->Draw("deltaPhiDegWithDynamicPitchRoll:zoomPhiDeg", "TMath::Abs(deltaPhiDegWithDynamicPitchRoll) < 3", "goff");
  // TH2F* hDeltaPhiDeg_vs_phiZoom2 = (TH2F*) htemp->Clone("hDeltaPhiDegWithDynamicPitchRoll_vs_phiZoom");
  // hDeltaPhiDeg_vs_phiZoom2->SetTitle("Close up of #phi_{expected} - #phi_{zoom} With dynamic correction; #phi_{zoom} (Degrees); #delta#phi (degrees); Events per bin");
  // hDeltaPhiDeg_vs_phiZoom2->Draw("colz");



  //   TCanvas* c3b = new TCanvas();
  //   chain->Draw("deltaThetaDegWithDynamicPitchRoll", "TMath::Abs(deltaThetaDegWithDynamicPitchRoll) < 3", "goff");
  //   TH1F* hDeltaThetaDeg2 = (TH1F*) htemp->Clone("hDeltaThetaDegWithDynamicPitchRoll");
  //   hDeltaThetaDeg2->SetTitle("Close up of #theta_{expected} - #theta_{zoom} With dynamic correction; #delta#theta (degrees); Events per bin");
  //   hDeltaThetaDeg2->Rebin(2);    
  //   hDeltaThetaDeg2->Draw();
  //   hDeltaThetaDeg2->GetYaxis()->SetNoExponent(1);  
  //   hDeltaThetaDeg2->Fit("gaus");


  
  //   TCanvas* c4b = new TCanvas();
  //   // chain->Draw("deltaThetaDeg:zoomPhiDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
  //   chain->Draw("deltaThetaDegWithDynamicPitchRoll:zoomPhiDeg", "TMath::Abs(deltaThetaDegWithDynamicPitchRoll) < 3", "goff"); 
  //   TH2F* hDeltaThetaDeg_vs_thetaZoom2 = (TH2F*) htemp->Clone("hDeltaThetaDegWithDynamicPitchRoll_vs_thetaZoom");
  //   hDeltaThetaDeg_vs_thetaZoom2->SetTitle("Close up of #theta_{expected} - #theta_{zoom} With dynamic correction; #phi_{zoom} (Degrees); #delta#theta (degrees); Events per bin");
  //   hDeltaThetaDeg_vs_thetaZoom2->Draw("colz");









    
    TCanvas* c1c = new TCanvas();
    chain->Draw("deltaPhiDegWithProvisionalBestConstantOffsets", "TMath::Abs(deltaPhiDegWithProvisionalBestConstantOffsets) < 3", "goff");
    TH1F* hDeltaPhiDeg3 = (TH1F*) htemp->Clone("hDeltaPhiDegWithProvisionalBestConstantOffsets");
    hDeltaPhiDeg3->SetTitle("Close up of #phi_{expected} - #phi_{zoom} Best fit offsets; #delta#phi (degrees); Events per bin");
    hDeltaPhiDeg3->Draw(); 
    hDeltaPhiDeg3->GetYaxis()->SetNoExponent(1);
    hDeltaPhiDeg3->Fit("gaus");

  
  
  TCanvas* c2c = new TCanvas();
  chain->Draw("deltaPhiDegWithProvisionalBestConstantOffsets:zoomPhiDeg", "TMath::Abs(deltaPhiDegWithProvisionalBestConstantOffsets) < 3", "goff");
  TH2F* hDeltaPhiDeg_vs_phiZoom3 = (TH2F*) htemp->Clone("hDeltaPhiDegWithProvisionalBestConstantOffsets_vs_phiZoom");
  hDeltaPhiDeg_vs_phiZoom3->SetTitle("Close up of #phi_{expected} - #phi_{zoom} Best fit offsets; #phi_{zoom} (Degrees); #delta#phi (degrees); Events per bin");
  hDeltaPhiDeg_vs_phiZoom3->Draw("colz");



    TCanvas* c3c = new TCanvas();
    chain->Draw("deltaThetaDegWithProvisionalBestConstantOffsets", "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 3", "goff");
    TH1F* hDeltaThetaDeg3 = (TH1F*) htemp->Clone("hDeltaThetaDegWithProvisionalBestConstantOffsets");
    hDeltaThetaDeg3->SetTitle("Close up of #theta_{expected} - #theta_{zoom} Best fit offsets; #delta#theta (degrees); Events per bin");
    hDeltaThetaDeg3->Rebin(2);    
    hDeltaThetaDeg3->Draw();
    hDeltaThetaDeg3->GetYaxis()->SetNoExponent(1);  
    hDeltaThetaDeg3->Fit("gaus");


  
    TCanvas* c4c = new TCanvas();
    // chain->Draw("deltaThetaDeg:zoomPhiDeg", "TMath::Abs(deltaThetaDeg) < 3", "goff");
    chain->Draw("deltaThetaDegWithProvisionalBestConstantOffsets:zoomPhiDeg", "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 3", "goff"); 
    TH2F* hDeltaThetaDeg_vs_thetaZoom3 = (TH2F*) htemp->Clone("hDeltaThetaDegWithProvisionalBestConstantOffsets_vs_thetaZoom");
    hDeltaThetaDeg_vs_thetaZoom3->SetTitle("Close up of #theta_{expected} - #theta_{zoom} Best fit offsets; #phi_{zoom} (Degrees); #delta#theta (degrees); Events per bin");
    hDeltaThetaDeg_vs_thetaZoom3->Draw("colz");
    

    
  }  
  
}
