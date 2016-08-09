#include "RootTools.h"
#include "ProgressBar.h"

void plotAngResVsSnr2(){

  TChain* cLdb = new TChain("angResTree");
  TChain* cSnrLdb = new TChain("snrTree");  

  TChain* cWais = new TChain("angResTree");
  TChain* cSnrWais = new TChain("snrTree");  
  
  // cLdb->Add("filterOffs/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
  // cSnrLdb->Add("filterOffs/phaseCenter/generateSignalToNoiseRatioTreeVPOLPlots_*");

  // cWais->Add("filterOffs/phaseCenter/generateAngularResolutionTreePlots_*");
  // cSnrWais->Add("filterOffs/phaseCenter/generateSignalToNoiseRatioTreePlots_*");

  cLdb->Add("extremeFilter/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
  cSnrLdb->Add("extremeFilter/phaseCenter/generateSignalToNoiseRatioTreeVPOLPlots_*");
  
  cWais->Add("filterOns/phaseCenter/generateAngularResolutionTreePlots_*");
  cSnrWais->Add("filterOns/phaseCenter/generateSignalToNoiseRatioTreePlots_*");
  
  Double_t deltaPhiDegLdb = 0;
  cLdb->SetBranchAddress("deltaPhiDeg", &deltaPhiDegLdb);
  Double_t deltaThetaDegLdb = 0;
  cLdb->SetBranchAddress("deltaThetaDeg", &deltaThetaDegLdb);
  UInt_t eventNumber = 0;
  cLdb->SetBranchAddress("eventNumber", &eventNumber);
  
  
  Double_t phiExpected = 0;
  cLdb->SetBranchAddress("phiExpected", &phiExpected);

  Double_t snr0Ldb = 0;
  cSnrLdb->SetBranchAddress("snr0", &snr0Ldb);

  Long64_t nEntries = cSnrLdb->GetEntries();
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

  const int nBinsPhi = 32;
  
  std::vector<double> theSnrBinsLdb = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21};//, 22, 23, 24, 25};

  TH2D* hSnrVsPhiExpectedLdb = new TH2D("hSnrVsPhiExpectedLdb", "SNR vs. #phi_{LDB}; #phi_{LDB} (Degrees); SNR (no units); Events per bin", nBinsPhi, 0, 360, theSnrBinsLdb.size()-1, &theSnrBinsLdb[0]);


  TH1D* hDeltaPhiLdb = new TH1D("hDeltaPhiLdb", "#delta#phi_{LDB}; #delta#phi (Degrees); Events per bin", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  // TH2D* hDeltaPhiLdb2 = new TH2D("hDeltaPhiLdb2", "#delta#phi_{LDB} vs. SNR; SNR (no units); #delta#phi (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiLdb2 = new TH2D("hDeltaPhiLdb2", "#delta#phi_{LDB} vs. SNR; SNR (no units); #delta#phi (Degrees); Events per bin", theSnrBinsLdb.size()-1, &theSnrBinsLdb[0], phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);  

  TH1D* hDeltaThetaLdb = new TH1D("hDeltaThetaLdb", "#delta#theta_{LDB}; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  // TH2D* hDeltaThetaLdb2 = new TH2D("hDeltaThetaLdb2", "#delta#theta_{LDB} vs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", nBinsSnr, minSnr, maxSnr, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaLdb2 = new TH2D("hDeltaThetaLdb2", "#delta#theta_{LDB} vs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", theSnrBinsLdb.size()-1, &theSnrBinsLdb[0], thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);  

  double lowestSnr = 9999;
  double highestSnr = -9999;
  UInt_t lowEventNumber;
  UInt_t highEventNumber;  
  
  const double phiCut = 3;  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cLdb->GetEntry(entry);
    cSnrLdb->GetEntry(entry);    
    // cLdbOff->GetEntry(entry);

    if(TMath::Abs(deltaPhiDegLdb) < phiCut){

      hSnrVsPhiExpectedLdb->Fill(phiExpected, snr0Ldb);
	
      // if(snr0Ldb > 6){
      hDeltaPhiLdb->Fill(deltaPhiDegLdb);
      hDeltaPhiLdb2->Fill(snr0Ldb, deltaPhiDegLdb);

      hDeltaThetaLdb->Fill(deltaThetaDegLdb);
      hDeltaThetaLdb2->Fill(snr0Ldb, deltaThetaDegLdb);
      // }

      if(snr0Ldb < lowestSnr){
	lowestSnr = snr0Ldb;
	lowEventNumber = eventNumber;
      }
      if(snr0Ldb > highestSnr){
	highestSnr = snr0Ldb;	
	highEventNumber = eventNumber;
      }
    }
    
    p.inc(entry, maxEntry);
  }

  std::cout << "LDB: " << std::endl;
  std::cout << lowEventNumber << "\t" << lowestSnr << std::endl;
  std::cout << highEventNumber << "\t" << highestSnr << std::endl;

  auto c0 = new TCanvas();
  hSnrVsPhiExpectedLdb->Draw("colz");
  
  TObjArray* tArrLdbPhi =  new TObjArray();
  hDeltaPhiLdb2->FitSlicesY(NULL, 1, hDeltaPhiLdb2->GetNbinsY(), 0, "QNR", tArrLdbPhi);
  TH1D* hResolutionVsSnrPhiLdb = NULL;
  for(int i=0; i < tArrLdbPhi->GetEntries(); i++){
    auto h = (TH1D*)tArrLdbPhi->At(i);
    // auto c0 = new TCanvas();
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrPhiLdb = h;
      hResolutionVsSnrPhiLdb->SetTitle("VPOL #delta#phi_{LDB}");      
    }
    new TCanvas();
    h->Draw();    
  }
  
  TObjArray* tArrLdbTheta =  new TObjArray();
  hDeltaThetaLdb2->FitSlicesY(NULL, 1, hDeltaThetaLdb2->GetNbinsY(), 0, "QNR", tArrLdbTheta);
  TH1D* hResolutionVsSnrThetaLdb = NULL;
  for(int i=0; i < tArrLdbTheta->GetEntries(); i++){
    auto h = (TH1D*)tArrLdbTheta->At(i);
    // auto c0 = new TCanvas();    
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrThetaLdb = h;
      hResolutionVsSnrThetaLdb->SetTitle("VPOL #delta#theta_{LDB}");
    }
    new TCanvas();
    h->Draw();
    
  }

  Double_t deltaPhiDegWais = 0;
  cWais->SetBranchAddress("deltaPhiDeg", &deltaPhiDegWais);
  Double_t deltaThetaDegWais = 0;
  cWais->SetBranchAddress("deltaThetaDeg", &deltaThetaDegWais);
  
  phiExpected = 0;
  cWais->SetBranchAddress("phiExpected", &phiExpected);

  Double_t snr0Wais = 0;
  cSnrWais->SetBranchAddress("snr0", &snr0Wais);

  std::vector<double> theSnrBinsWais = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 18, 21};//, 22, 23, 24, 25};


  TH2D* hSnrVsPhiExpectedWais = new TH2D("hSnrVsPhiExpectedWais", "SNR vs. #phi_{WAIS}; #phi_{WAIS} (Degrees); SNR (no units); Events per bin", nBinsPhi, 0, 360, theSnrBinsWais.size()-1, &theSnrBinsWais[0]);

    
  TH1D* hDeltaPhiWais = new TH1D("hDeltaPhiWais", "#delta#phi_{WAIS}; #delta#phi (Degrees); Events per bin", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiWais2 = new TH2D("hDeltaPhiWais2", "#delta#phi_{WAIS} vs. SNR; SNR (no units); #delta#phi (Degrees); Events per bin", theSnrBinsWais.size()-1, &theSnrBinsWais[0], phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH1D* hDeltaThetaWais = new TH1D("hDeltaThetaWais", "#delta#theta_{WAIS}; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaWais2 = new TH2D("hDeltaThetaWais2", "#delta#theta_{WAIS} vs. SNR; SNR (no units); #delta#theta (Degrees); Events per bin", theSnrBinsWais.size()-1, &theSnrBinsWais[0], thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  Long64_t nEntries2 = cSnrWais->GetEntries();
  Long64_t maxEntry2 = 0; //5; //1000; //10000;
  Long64_t startEntry2 = 0;
  if(maxEntry2<=0 || maxEntry2 > nEntries2) maxEntry2 = nEntries2;
  std::cout << "Processing " << maxEntry2 << " of " << nEntries2 << " entries." << std::endl;
  ProgressBar p2(maxEntry2-startEntry2);
  
  for(Long64_t entry = startEntry2; entry < maxEntry2; entry++){
    cWais->GetEntry(entry);
    cSnrWais->GetEntry(entry);    

    if(TMath::Abs(deltaPhiDegWais) < phiCut){
      hSnrVsPhiExpectedWais->Fill(phiExpected, snr0Wais);
      hDeltaPhiWais->Fill(deltaPhiDegWais);
      hDeltaPhiWais2->Fill(snr0Wais, deltaPhiDegWais);

      hDeltaThetaWais->Fill(deltaThetaDegWais);
      hDeltaThetaWais2->Fill(snr0Wais, deltaThetaDegWais);
    }
    
    p2.inc(entry, maxEntry2);
  }  


  TObjArray* tArrWaisPhi =  new TObjArray();
  hDeltaPhiWais2->FitSlicesY(NULL, 1, hDeltaPhiWais2->GetNbinsY(), 0, "QNR", tArrWaisPhi);
  TH1D* hResolutionVsSnrPhiWais = NULL;
  for(int i=0; i < tArrWaisPhi->GetEntries(); i++){
    auto h = (TH1D*)tArrWaisPhi->At(i);
    // auto c0 = new TCanvas();
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrPhiWais = h;
      hResolutionVsSnrPhiWais->SetTitle("HPOL #delta#phi_{WAIS}");
    }
    new TCanvas();
    h->Draw();
    
  }
  
  TObjArray* tArrWaisTheta =  new TObjArray();
  hDeltaThetaWais2->FitSlicesY(NULL, 1, hDeltaThetaWais2->GetNbinsY(), 0, "QNR", tArrWaisTheta);
  TH1D* hResolutionVsSnrThetaWais = NULL;
  for(int i=0; i < tArrWaisTheta->GetEntries(); i++){
    auto h = (TH1D*)tArrWaisTheta->At(i);
    // auto c0 = new TCanvas();    
    // TString opt = i==0? "" : "same";
    // h->Draw(opt);
    if(i==2){
      hResolutionVsSnrThetaWais = h;
      hResolutionVsSnrThetaWais->SetTitle("HPOL #delta#theta_{WAIS}");      
    }
    new TCanvas();
    h->Draw();
  }
  
  
  auto c0a = new TCanvas();
  hSnrVsPhiExpectedWais->Draw("colz");
  return;
  
  auto c1 = new TCanvas();
  hDeltaPhiWais2->Draw("colz");

  auto c1a = new TCanvas();
  hDeltaPhiLdb2->Draw("colz");

  auto c2 = new TCanvas();
  hDeltaThetaWais2->Draw("colz");
  
  auto c2a = new TCanvas();
  hDeltaThetaLdb2->Draw("colz");

  // auto c3 = new TCanvas();
  // // hResolutionVsSnrPhiLdb->SetLineColor(kRed);
  // // hResolutionVsSnrThetaLdb->SetLineColor(kBlue);
  // hResolutionVsSnrPhiLdb->Draw();
  // hResolutionVsSnrThetaLdb->Draw("same");

  // auto c3a = new TCanvas();
  // // hResolutionVsSnrPhiWais->SetLineColor(kRed);
  // // hResolutionVsSnrThetaWais->SetLineColor(kBlue);
  // hResolutionVsSnrPhiWais->Draw();
  // hResolutionVsSnrThetaWais->Draw("same");

  auto c4 = new TCanvas();
  hResolutionVsSnrPhiWais->SetLineColor(kRed);
  hResolutionVsSnrPhiLdb->SetLineColor(kBlue);
  hResolutionVsSnrThetaWais->SetLineColor(kMagenta);
  hResolutionVsSnrThetaLdb->SetLineColor(kCyan);

  hResolutionVsSnrPhiWais->SetMarkerSize(0);
  hResolutionVsSnrPhiLdb->SetMarkerSize(0);
  hResolutionVsSnrThetaWais->SetMarkerSize(0);
  hResolutionVsSnrThetaLdb->SetMarkerSize(0);
  
  hResolutionVsSnrPhiLdb->Draw();
  hResolutionVsSnrThetaLdb->Draw("same");  
  hResolutionVsSnrPhiWais->Draw("same");
  hResolutionVsSnrThetaWais->Draw("same");

  
  auto l4 = c4->BuildLegend();
  // l1->SetNColumns(2);
  l4->Draw();
  // hs[0]->SetMaximum(0.13);

  hResolutionVsSnrPhiLdb->SetTitle("Angular Resolution vs. SNR; SNR; Angular Resolution (Degrees)");
  hResolutionVsSnrPhiLdb->SetMaximum(1);
  hResolutionVsSnrPhiLdb->SetMinimum(0);  
  
  //hDeltaPhiDegPhaseLdb->SetTitle("ANITA Angular Resolution Using Fitted Phase Centres; #delta#theta or #delta#phi (Degrees); Fraction of events per bin");

  
}
