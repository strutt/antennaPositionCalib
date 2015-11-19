void drawAngularResolutionPlots4(){

  gStyle->SetOptFit(01111);
  
  Int_t colCorrected = 2;
  Int_t colNotCorrected = 4;  
  Bool_t friends = kTRUE;
  auto chain = new TChain("angResTree");
  auto chain2 = new TChain("angResTreeFriend");

  for(int run=331; run<=354; run++){
    // TString fileName = TString::Format("photogrammetryNoPitchRoll/generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    // TString fileName = TString::Format("generateAngularResolutionTreeFriend_run%d-%dPlots.root", run, run);
    TString fileName = TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run);    
    chain->Add(fileName);
    // fileName = TString::Format("generateAngularResolutionTree_run%d-%dLindaPlots.root", run, run);
    // // fileName = TString::Format("generateAngularResolutionTreeFriend_run%d-%dPlots.root", run, run);
    // chain2->Add(fileName);
  }
  // chain->AddFriend(chain2);
  // chain->Show(0);
  // return;

  auto c5 = new TCanvas();
  chain->Draw("deltaPhiDeg0:phiExpected0", "TMath::Abs(deltaPhiDeg0) < 3", "goff");
  TH2F* hDeltaPhiDeg0 = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDeg0");
  hDeltaPhiDeg0->Draw("colz");

  auto c6 = new TCanvas();
  chain->Draw("deltaThetaDeg0:phiExpected0", "TMath::Abs(deltaThetaDeg0) < 3", "goff");
  TH2F* hDeltaThetaDeg0 = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDeg0");
  hDeltaThetaDeg0->Draw("colz");

  auto c7 = new TCanvas();
  chain->Draw("deltaPhiDegTilted:phiExpectedTilted", "TMath::Abs(deltaPhiDegTilted) < 3", "goff");
  TH2F* hDeltaPhiDegTilted = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDegTilted");
  hDeltaPhiDegTilted->Draw("colz");

  auto c8 = new TCanvas();
  chain->Draw("deltaThetaDegTilted:phiExpectedTilted", "TMath::Abs(deltaThetaDegTilted) < 3", "goff");
  TH2F* hDeltaThetaDegTilted = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDegTilted");
  hDeltaThetaDegTilted->Draw("colz");
  

  // auto c1 = new TCanvas();
  // chain->Draw("deltaPhiDegTilted:zoomPhiDeg", "TMath::Abs(deltaPhiDegTilted) < 3", "goff");
  // TH2F* hDeltaPhiDegTilted2D = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDegTilted2D");
  // hDeltaPhiDegTilted2D->SetTitle("Close up of #phi_{expected} - #phi_{zoom} Best heading/pitch/roll offset; Azimuth (degrees); #delta#phi (degrees); Events per bin");
  // hDeltaPhiDegTilted2D->Draw("colz");

  // auto c3 = new TCanvas();
  // chain->Draw("deltaPhiDeg0:zoomPhiDeg", "TMath::Abs(deltaPhiDeg0) < 3", "goff");
  // TH2F* hDeltaPhiDeg0_2D = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaPhiDeg0_2D");
  // hDeltaPhiDeg0_2D->SetTitle("Close up of #phi_{expected} - #phi_{zoom} No pitch roll offsets; Azimuth (degrees); #delta#phi (degrees); Events per bin");
  // hDeltaPhiDeg0_2D->Draw("colz");

  // auto c2 = new TCanvas();
  // chain->Draw("deltaThetaDegTilted:zoomPhiDeg", "TMath::Abs(deltaThetaDegTilted) < 3", "goff");
  // TH2F* hDeltaThetaDegTilted2D = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDegTilted2D");
  // hDeltaThetaDegTilted2D->SetTitle("Close up of #theta_{expected} - #theta_{zoom} Best heading/pitch/roll offset; Azimuth (degrees); #delta#theta (degrees); Events per bin");
  // hDeltaThetaDegTilted2D->Draw("colz");
  // hDeltaThetaDegTilted2D->GetYaxis()->SetNoExponent(1);

  // auto c4 = new TCanvas();
  // chain->Draw("deltaThetaDeg0:zoomPhiDeg", "TMath::Abs(deltaThetaDeg0) < 3", "goff");
  // TH2F* hDeltaThetaDeg0_2D = (TH2F*) gROOT->Get("htemp")->Clone("hDeltaThetaDeg0_2D");
  // hDeltaThetaDeg0_2D->SetTitle("Close up of #theta_{expected} - #theta_{zoom} No pitch roll offsets; Azimuth (degrees); #delta#theta (degrees); Events per bin");
  // hDeltaThetaDeg0_2D->Draw("colz");
  // hDeltaThetaDeg0_2D->GetYaxis()->SetNoExponent(1);

  
}
