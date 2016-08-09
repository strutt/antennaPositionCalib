#include "RootTools.h"
#include "ProgressBar.h"

void plotAngResVsSnr(){

  TChain* cFilt = new TChain("angResTree");
  // TChain* cFiltOff = new TChain("angResTree");
  TChain* cSnr = new TChain("snrTree");  

  // cFilt->Add("filterOns/phaseCenter/generateAngularResolutionTreePlots_*");
  // cFilt->Add("filterOffs/phaseCenter/generateAngularResolutionTreePlots_*");
  // cSnr->Add("filterOffs/phaseCenter/generateSignalToNoiseRatioTreePlots_*");

  // cFilt->Add("filterOns/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
  // cSnr->Add("filterOns/phaseCenter/generateSignalToNoiseRatioTreeVPOLPlots_*");
  cFilt->Add("filterOffs/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");  
  cSnr->Add("filterOffs/phaseCenter/generateSignalToNoiseRatioTreeVPOLPlots_*");

    
  Double_t zoomPhiDeg = 0;
  cFilt->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
  Double_t deltaPhiDeg = 0;  
  cFilt->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);  
  // Double_t zoomPhiDegOff = 0;
  // cFiltOff->SetBranchAddress("zoomPhiDeg", &zoomPhiDegOff);
  Double_t zoomThetaDeg = 0;
  cFilt->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);
  Double_t deltaThetaDeg = 0;  
  cFilt->SetBranchAddress("deltaThetaDeg", &deltaThetaDeg);  
  // Double_t zoomThetaDegOff = 0;
  // cFiltOff->SetBranchAddress("zoomThetaDeg", &zoomThetaDegOff);
  
  Double_t phiExpected = 0;
  cFilt->SetBranchAddress("phiExpected", &phiExpected);
  Double_t heading = 0;
  cFilt->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  cFilt->SetBranchAddress("realTime", &realTime);

  Double_t snr0 = 0;
  cSnr->SetBranchAddress("snr0", &snr0);

  Long64_t nEntries = cSnr->GetEntries();
  Long64_t maxEntry = 0; //5; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  double dPhi = 0.05;
  double phiRange = 6;

  double dTheta = 0.05;
  double thetaRange = 6;
  double rebin = 2;

  const int nBinsSnr = 64; //32
  // const int nBinsSnr = 32;
  const double minSnr = 0;
  const double maxSnr = 20; //10

  TH1D* hPhiDiff = new TH1D("hPhiDiff", "#phi_{on} - #phi_{off}; #delta#phi (Degrees); Events per bin", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hPhiDiff2 = new TH2D("hPhiDiff2", "#phi_{on} - #phi_{off}; SNR (no units); #delta#phi (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH1D* hThetaDiff = new TH1D("hThetaDiff", "#delta#theta vs. SNR; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hThetaDiff2 = new TH2D("hThetaDiff2", "#delta#theta cs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  const double phiCut = 3;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cFilt->GetEntry(entry);
    cSnr->GetEntry(entry);    
    // cFiltOff->GetEntry(entry);

    if(TMath::Abs(deltaPhiDeg) < phiCut){

      hPhiDiff->Fill(deltaPhiDeg);
      hPhiDiff2->Fill(snr0, deltaPhiDeg);

      hThetaDiff->Fill(deltaThetaDeg);
      hThetaDiff2->Fill(snr0, deltaThetaDeg);
    }
    
    p.inc(entry, maxEntry);
  }

  TObjArray* tArrPhi =  new TObjArray();
  hPhiDiff2->FitSlicesY(NULL, 1, hPhiDiff2->GetNbinsY(), 0, "QNR", tArrPhi);
  TH1D* hResolutionVsSnrPhi = NULL;
  for(int i=0; i < tArrPhi->GetEntries(); i++){
    auto h = (TH1D*)tArrPhi->At(i);
    // auto c0 = new TCanvas();
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrPhi = h;
    }
  }
  
  TObjArray* tArrTheta =  new TObjArray();
  hThetaDiff2->FitSlicesY(NULL, 1, hThetaDiff2->GetNbinsY(), 0, "QNR", tArrTheta);
  TH1D* hResolutionVsSnrTheta = NULL;
  for(int i=0; i < tArrTheta->GetEntries(); i++){
    auto h = (TH1D*)tArrTheta->At(i);
    // auto c0 = new TCanvas();    
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrTheta = h;
    }
  }
  
  // auto c1 = new TCanvas();
  // // c1->SetLogy(1);
  // hPhiDiff->Draw();

  auto c1a = new TCanvas();
  // c1->SetLogy(1);
  hPhiDiff2->Draw("colz");

  // auto c2 = new TCanvas();
  // // c1->SetLogy(1);
  // hThetaDiff->Draw();
  
  auto c2a = new TCanvas();
  // c1->SetLogy(1);
  hThetaDiff2->Draw("colz");

  auto c3 = new TCanvas();
  hResolutionVsSnrPhi->SetLineColor(kRed);
  hResolutionVsSnrTheta->SetLineColor(kBlue);
  hResolutionVsSnrPhi->Draw();
  hResolutionVsSnrTheta->Draw("same");

  
}
