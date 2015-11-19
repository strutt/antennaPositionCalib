void drawAngularResolutionPlots3(){

  gStyle->SetOptFit(01111);
  
  auto chain = new TChain("angResTree");
  auto chain2 = new TChain("angResTreeWithDynamicPitchRoll");  

  chain->Add("generateAngularResolutionTree_run*Plots.root");
  chain2->Add("generateAngularResolutionTreeFriend_run*Plots.root");
  chain->AddFriend(chain2);

  // TString cuts = "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 2 && TMath::Abs(roll+0.5) < 1";
  TString cuts = "TMath::Abs(deltaThetaDegWithProvisionalBestConstantOffsets) < 2 && TMath::Abs(roll+0.5) < 1.5 && TMath::Abs(pitch-0.5) < 2";  
  auto c0 = new TCanvas();

  chain->Draw("pitch>>hPitch", cuts, "goff");
  auto hPitch = (TH1F*) gROOT->Get("hPitch");
  hPitch->SetLineColor(kRed);
  hPitch->Rebin(2);  
  hPitch->SetLineStyle(2);  
  hPitch->Draw();

  chain->Draw("roll>>hRoll", cuts, "goff");
  auto hRoll = (TH1F*) gROOT->Get("hRoll");
  hRoll->SetLineColor(kBlue);
  hRoll->SetLineStyle(2);
  hRoll->Rebin(2);    
  hRoll->Draw("same");

  cuts += " && (!(roll==0 && pitch==0)) ";

  chain->Draw("pitch>>hPitch2", cuts, "goff");
  auto hPitch2 = (TH1F*) gROOT->Get("hPitch2");
  hPitch2->SetLineColor(kRed);
  hPitch2->Draw("same");
  hPitch2->Rebin(2);
  
  chain->Draw("roll>>hRoll2", cuts, "goff");
  auto hRoll2 = (TH1F*) gROOT->Get("hRoll2");
  hRoll2->SetLineColor(kBlue);
  hRoll2->Rebin(2);
  hRoll2->Draw("same");
  // c0->SetLogy(1);

  cout << "pitch rms = " << hPitch2->GetRMS() << endl;  
  cout << "roll rms = " << hRoll2->GetRMS() << endl;

  hPitch->SetTitle("Dynamic pitch and roll information for WAIS pulser events");
  hPitch->GetXaxis()->SetTitle("Dynamic offset (degrees)");
  hPitch->GetYaxis()->SetTitle("Events per bin");
  auto l0 = new TLegend(0.8, 0.8, 0.98, 1);
  l0->AddEntry(hPitch2, "Pitch", "l");
  l0->AddEntry(hRoll2, "Roll", "l");
  l0->Draw();
  
  auto c1 = new TCanvas();
  chain->Draw("pitch:roll>>hPitchRoll", cuts, "colz");
  auto hPitchRoll = (TH2D*) gROOT->Get("hPitchRoll");
  auto hPitchRoll_pfx = hPitchRoll->ProfileX();
  c1->SetLogz(1);
  hPitchRoll_pfx->SetLineWidth(2);
  hPitchRoll_pfx->Draw("same");
  // TF1* fPitchRoll = new TF1("fPitchRoll", "[0]*x*x + [1]*x +[2]");
  // fPitchRoll->SetLineWidth(2);
  // fPitchRoll->SetLineColor(kMagenta);
  // hPitchRoll_pfx->Fit(fPitchRoll);

  // auto l1 = new TLegend(0.8, 0.8, 0.98, 1);
  // l1->AddEntry(hPitchRoll_pfx, "Profile", "l");
  // // l1->AddEntry(fPitchRoll, "Linear fit to profile", "l");
  // l1->Draw();

  hPitchRoll->SetTitle("Dynamic Pitch vs. Dynamic Roll");
  hPitchRoll->GetXaxis()->SetTitle("Roll (Degrees)");
  hPitchRoll->GetYaxis()->SetTitle("Pitch (Degrees)");
  cout << hPitchRoll->GetCorrelationFactor() << endl;
  c1->SetLogz(1);
  


  auto c1_n2 = new TCanvas();
  chain->Draw("roll:deltaThetaDegWithProvisionalBestConstantOffsets>>hCorrectedRoll", cuts, "colz");
  auto hCorrectedRoll = (TH2D*) gROOT->Get("hCorrectedRoll");
  auto hCorrectedRoll_pfx = hCorrectedRoll->ProfileX();
  c1_n2->SetLogz(1);
  hCorrectedRoll_pfx->SetLineWidth(2);
  hCorrectedRoll_pfx->Draw("same");
  // TF1* fCorrectedRoll = new TF1("fCorrectedRoll", "[0]*x*x + [1]*x +[2]");
  // fCorrectedRoll->SetLineWidth(2);
  // fCorrectedRoll->SetLineColor(kMagenta);
  // // hCorrectedRoll->Fit(fCorrectedRoll);
  hCorrectedRoll->SetTitle("#delta#theta vs Roll (static best fit offsets applied)");
  hCorrectedRoll->GetXaxis()->SetTitle("#delta#theta (Degrees)");
  hCorrectedRoll->GetYaxis()->SetTitle("Roll (Degrees)");
  cout << hCorrectedRoll->GetCorrelationFactor() << endl;


  auto c1_n3 = new TCanvas();
  chain->Draw("pitch:deltaThetaDegWithProvisionalBestConstantOffsets>>hCorrectedPitch", cuts, "colz");
  auto hCorrectedPitch = (TH2D*) gROOT->Get("hCorrectedPitch");
  auto hCorrectedPitch_pfx = hCorrectedPitch->ProfileX();
  c1_n3->SetLogz(1);
  hCorrectedPitch_pfx->SetLineWidth(2);
  hCorrectedPitch_pfx->Draw("same");
  // TF1* fCorrectedPitch = new TF1("fCorrectedPitch", "[0]*x*x + [1]*x +[2]");
  // fCorrectedPitch->SetLineWidth(2);
  // fCorrectedPitch->SetLineColor(kMagenta);
  // hCorrectedPitch->Fit(fCorrectedPitch);
  hCorrectedPitch->SetTitle("#delta#theta vs Pitch (static best fit offsets applied)");
  hCorrectedPitch->GetXaxis()->SetTitle("#delta#theta (Degrees)");
  hCorrectedPitch->GetYaxis()->SetTitle("Pitch (Degrees)");
  cout << hCorrectedPitch->GetCorrelationFactor() << endl;

  
  return;
  
  auto c2 = new TCanvas();
  chain->Draw("pitch:deltaThetaDeg", cuts, "colz");  
  
  auto c3 = new TCanvas();
  chain->Draw("roll:deltaThetaDeg", cuts, "colz");

  auto c4 = new TCanvas();
  chain->Draw("pitch+roll:deltaThetaDeg", cuts, "colz");


  auto c5 = new TCanvas();
  chain->Draw("deltaThetaDeg:pitch+roll", "TMath::Abs(deltaThetaDeg) < 2 && TMath::Abs(roll+0.5) < 1 && roll!=0 && zoomPhiDeg >= 0 && zoomPhiDeg < 360", "colz");

  gROOT->ProcessLine("TProfile* hprof = htemp->ProfileX(); hprof->Draw(\"same\"); hprof->SetLineWidth(2)");
  

  
  // chain->Draw("pitch+roll", cuts, "same"); //"colz");

  
  // auto c5 = new TCanvas();
  // chain->Draw("deltaThetaDeg:zoomPhiDeg", "TMath::Abs(deltaThetaDeg) < 2 && TMath::Abs(roll+0.5) < 1 && roll!=0 && zoomPhiDeg < 90", "colz");
  
}
