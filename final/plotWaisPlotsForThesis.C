#include "RootTools.h"
#include "ProgressBar.h"

void plotWaisPlotsForThesis(){

  TChain* cPhoto = new TChain("angResTree");
  TChain* cPhase = new TChain("angResTree");
  TChain* cSnrPhoto = new TChain("snrTree");
  TChain* cSnrPhase = new TChain("snrTree");

  cPhase->Add("filterOffs/phaseCenter/generateAngularResolutionTreePlots_*");
  cPhoto->Add("filterOffs/photogrammetry/generateAngularResolutionTreePlots_*");
  cSnrPhase->Add("filterOffs/phaseCenter/generateSignalToNoiseRatioTreePlots_*");
  cSnrPhoto->Add("filterOffs/photogrammetry/generateSignalToNoiseRatioTreePlots_*");
    
  Double_t deltaPhiDegPhase = 0;  
  cPhase->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhase);
  Double_t deltaPhiDegPhoto = 0;  
  cPhoto->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhoto);

  Double_t deltaThetaDegPhase = 0;  
  cPhase->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhase);
  Double_t deltaThetaDegPhoto = 0;  
  cPhoto->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhoto);  
  
  
  
  Double_t phiExpected = 0;
  cPhase->SetBranchAddress("phiExpected", &phiExpected);
  Double_t heading = 0;
  cPhase->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  cPhase->SetBranchAddress("realTime", &realTime);

  Double_t snr0Phase = 0;
  cSnrPhase->SetBranchAddress("snr0", &snr0Phase);
  Double_t snr0Photo = 0;
  cSnrPhoto->SetBranchAddress("snr0", &snr0Photo);

  Long64_t nEntries = cSnrPhase->GetEntries();
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

  const int nBinsPhi = 256;

  
  TH1D* hDeltaPhiDegPhase = new TH1D("hDeltaPhiDegPhase", "HPol #phi_{on} - #phi_{off}; #delta#phi (Degrees); Events per bin", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhase2 = new TH2D("hDeltaPhiDegPhase2", "HPol #phi_{on} - #phi_{off}; SNR (no units); #delta#phi (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhase3 = new TH2D("hDeltaPhiDegPhase3", "HPol #phi_{WAIS} - #phi_{measured} Fitted Phase Centre Antenna Positions; #phi_{WAIS} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH1D* hDeltaThetaDegPhase = new TH1D("hDeltaThetaDegPhase", "HPol #delta#theta vs. SNR; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhase2 = new TH2D("hDeltaThetaDegPhase2", "HPol #delta#theta cs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhase3 = new TH2D("hDeltaThetaDegPhase3", "HPol #theta_{WAIS} - #theta_{measured} Fitted Phase Centre Antenna Positions; #theta_{WAIS} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  
  TH1D* hDeltaPhiDegPhoto = new TH1D("hDeltaPhiDegPhoto", "HPol #phi_{on} - #phi_{off}; #delta#phi (Degrees); Events per bin", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhoto2 = new TH2D("hDeltaPhiDegPhoto2", "HPol #phi_{on} - #phi_{off}; SNR (no units); #delta#phi (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhoto3 = new TH2D("hDeltaPhiDegPhoto3", "HPol #phi_{WAIS} - #phi_{measured} Photogrammetry Antenna Positions; #phi_{WAIS} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH1D* hDeltaThetaDegPhoto = new TH1D("hDeltaThetaDegPhoto", "HPol #delta#theta vs. SNR; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhoto2 = new TH2D("hDeltaThetaDegPhoto2", "HPol #delta#theta cs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhoto3 = new TH2D("hDeltaThetaDegPhoto3", "HPol #theta_{WAIS} - #theta_{measured} Photogrammetry Antenna Positions; #theta_{WAIS} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  
  const double phiCut = 3;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cPhase->GetEntry(entry);
    cPhoto->GetEntry(entry);    
    cSnrPhase->GetEntry(entry);
    cSnrPhoto->GetEntry(entry);        
    // cFiltOff->GetEntry(entry);

    if(TMath::Abs(deltaPhiDegPhase) < phiCut){

      hDeltaPhiDegPhase->Fill(deltaPhiDegPhase);
      hDeltaThetaDegPhase->Fill(deltaThetaDegPhase);
      hDeltaPhiDegPhase2->Fill(snr0Phase, deltaPhiDegPhase);
      hDeltaThetaDegPhase2->Fill(snr0Phase, deltaThetaDegPhase);
      hDeltaPhiDegPhase3->Fill(phiExpected, deltaPhiDegPhase);
      hDeltaThetaDegPhase3->Fill(phiExpected, deltaThetaDegPhase);

      
    }
    if(TMath::Abs(deltaPhiDegPhoto) < phiCut){
      hDeltaPhiDegPhoto->Fill(deltaPhiDegPhoto);
      hDeltaThetaDegPhoto->Fill(deltaThetaDegPhoto);
      hDeltaPhiDegPhoto2->Fill(snr0Photo, deltaPhiDegPhoto);
      hDeltaThetaDegPhoto2->Fill(snr0Photo, deltaThetaDegPhoto);
      hDeltaPhiDegPhoto3->Fill(phiExpected, deltaPhiDegPhoto);
      hDeltaThetaDegPhoto3->Fill(phiExpected, deltaThetaDegPhoto);
    }
    
    p.inc(entry, maxEntry);
  }

  TObjArray* tArrPhi =  new TObjArray();
  hDeltaPhiDegPhase2->FitSlicesY(NULL, 1, hDeltaPhiDegPhase2->GetNbinsY(), 0, "QNR", tArrPhi);
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
  hDeltaThetaDegPhase2->FitSlicesY(NULL, 1, hDeltaThetaDegPhase2->GetNbinsY(), 0, "QNR", tArrTheta);
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

  auto c0 = new TCanvas();
  hDeltaPhiDegPhoto3->Draw("colz");
  auto c0a = new TCanvas();
  hDeltaThetaDegPhoto3->Draw("colz");
  
  
  auto c0b = new TCanvas();
  hDeltaPhiDegPhase3->Draw("colz");
  auto c0c = new TCanvas();
  hDeltaThetaDegPhase3->Draw("colz");
  
  
  
  auto c1 = new TCanvas();
  // c1->SetLogy(1);
  hDeltaPhiDegPhase->SetLineColor(kRed);
  hDeltaPhiDegPhoto->SetLineColor(kBlue);
  hDeltaPhiDegPhase->Draw();  
  hDeltaPhiDegPhoto->Draw("same");

  auto c1a = new TCanvas();
  // c1->SetLogy(1);
  hDeltaPhiDegPhase2->Draw("colz");

  auto c2 = new TCanvas();
  // c1->SetLogy(1);
  hDeltaThetaDegPhase->SetLineColor(kRed);
  hDeltaThetaDegPhoto->SetLineColor(kBlue);
  hDeltaThetaDegPhase->Draw();  
  hDeltaThetaDegPhoto->Draw("same");
  
  auto c2a = new TCanvas();
  // c1->SetLogy(1);
  hDeltaThetaDegPhase2->Draw("colz");

  auto c3 = new TCanvas();
  hResolutionVsSnrPhi->SetLineColor(kRed);
  hResolutionVsSnrTheta->SetLineColor(kBlue);
  hResolutionVsSnrPhi->Draw();
  hResolutionVsSnrTheta->Draw("same");

}
