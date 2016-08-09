#include "RootTools.h"
#include "ProgressBar.h"

void drawAngRes(){

  TChain* c = new TChain("angResTree");
  c->Add("extremeFilter/phaseCenter/generateAngularResolutionTreeVPOLPlots_1*");

  Double_t deltaPhiDeg = 0;
  c->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);


  Double_t deltaThetaDeg = 0;
  c->SetBranchAddress("deltaThetaDeg", &deltaThetaDeg);
  
  
  Double_t phiExpected = 0;
  c->SetBranchAddress("phiExpected", &phiExpected);
  Double_t heading = 0;
  c->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  c->SetBranchAddress("realTime", &realTime);

  
  Long64_t nEntries = c->GetEntries();
  Long64_t maxEntry = 0; //5; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  double dPhi = 0.05;
  double phiRange = 6;

  double dTheta = 0.05;
  double thetaRange = 6;


  
  TH1D* hDeltaPhiDeg = new TH1D("hDeltaPhiDeg", "#phi_{on} - #phi_{off}; #delta#phi (Degrees); Events per bin", phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH1D* hDeltaThetaDeg = new TH1D("hDeltaThetaDeg", "#theta_{on} - #theta_{off}; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  


  const int nBinsPhi = 360;
  TH2D* hDeltaPhiDeg2 = new TH2D("hDeltaPhiDeg2", "#phi_{on} - #phi_{off}; #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH2D* hDeltaThetaDeg2 = new TH2D("hDeltaThetaDeg2", "#theta_{on} - #theta_{off}; #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    c->GetEntry(entry);

    hDeltaPhiDeg->Fill(deltaPhiDeg);
    hDeltaThetaDeg->Fill(deltaThetaDeg);
    hDeltaPhiDeg2->Fill(phiExpected, deltaPhiDeg);
    hDeltaThetaDeg2->Fill(phiExpected, deltaThetaDeg);
    
    p.inc(entry, maxEntry);
  }
  auto c1 = new TCanvas();
  hDeltaPhiDeg->Draw();

  auto c2 = new TCanvas();
  hDeltaPhiDeg2->Draw("colz");

  auto c3 = new TCanvas();
  hDeltaThetaDeg->Draw();

  auto c4 = new TCanvas();
  hDeltaThetaDeg2->Draw("colz");
  
}
