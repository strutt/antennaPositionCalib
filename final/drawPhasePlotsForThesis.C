#include "RootTools.h"
#include "ProgressBar.h"

void plotPhotoPlotsForThesis(){

  TChain* cPhotoLdb = new TChain("angResTree");
  TChain* cPhotoWais = new TChain("angResTree");

  cPhotoWais->Add("filterOffs/photogrammetry/generateAngularResolutionTreePlots_*");
  cPhotoLdb->Add("filterOffs/photogrammetry/generateAngularResolutionTreeVPOLPlots_*");
    
  Double_t deltaPhiDegPhotoWais = 0;  
  cPhotoWais->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhotoWais);
  Double_t deltaPhiDegPhotoLdb = 0;  
  cPhotoLdb->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhotoLdb);

  Double_t deltaThetaDegPhotoWais = 0;  
  cPhotoWais->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhotoWais);
  Double_t deltaThetaDegPhotoLdb = 0;  
  cPhotoLdb->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhotoLdb);  
  
  
  
  Double_t phiExpected = 0;
  cPhotoWais->SetBranchAddress("phiExpected", &phiExpected);
  Double_t heading = 0;
  cPhotoWais->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  cPhotoWais->SetBranchAddress("realTime", &realTime);

  Long64_t nEntries = cPhotoWais->GetEntries();
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

  
  TH1D* hDeltaPhiDegPhotoWais = new TH1D("hDeltaPhiDegPhotoWais", "HPol #delta#phi_{WAIS}", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhotoWais3 = new TH2D("hDeltaPhiDegPhotoWais3", "HPol #delta#phi_{WAIS}; #delta#phi_{WAIS} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH1D* hDeltaThetaDegPhotoWais = new TH1D("hDeltaThetaDegPhotoWais", "HPol #delta#theta_{WAIS}", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhotoWais3 = new TH2D("hDeltaThetaDegPhotoWais3", "HPol #delta#theta_{WAIS} Fitted PhotoWais Centre Antenna Positions; #delta#theta_{WAIS} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);


  
  TH1D* hDeltaPhiDegPhotoLdb = new TH1D("hDeltaPhiDegPhotoLdb", "VPol #delta#phi_{LDB}", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);  
  TH2D* hDeltaPhiDegPhotoLdb3 = new TH2D("hDeltaPhiDegPhotoLdb3", "HPol #delta#phi_{WAIS} PhotoLdbgrammetry Antenna Positions; #delta#phi_{WAIS} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH1D* hDeltaThetaDegPhotoLdb = new TH1D("hDeltaThetaDegPhotoLdb", "VPol #delta#theta_{LDB}", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);  
  TH2D* hDeltaThetaDegPhotoLdb3 = new TH2D("hDeltaThetaDegPhotoLdb3", "HPol #delta#theta_{WAIS} PhotoLdbgrammetry Antenna Positions; #delta#theta_{WAIS} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  
  const double phiCut = 3;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cPhotoWais->GetEntry(entry);

    if(TMath::Abs(deltaPhiDegPhotoWais) < phiCut){

      hDeltaPhiDegPhotoWais->Fill(deltaPhiDegPhotoWais);
      hDeltaPhiDegPhotoWais3->Fill(phiExpected, deltaPhiDegPhotoWais);

      hDeltaThetaDegPhotoWais->Fill(deltaThetaDegPhotoWais);      
      hDeltaThetaDegPhotoWais3->Fill(phiExpected, deltaThetaDegPhotoWais);
      
    }
    p.inc(entry, maxEntry);    
  }

  Long64_t nEntries2 = cPhotoLdb->GetEntries();
  Long64_t maxEntry2 = 0; //5; //1000; //10000;
  Long64_t startEntry2 = 0;
  if(maxEntry2<=0 || maxEntry2 > nEntries2) maxEntry2 = nEntries2;
  std::cout << "Processing " << maxEntry2 << " of " << nEntries2 << " entries." << std::endl;
  ProgressBar p2(maxEntry2-startEntry2);
  
  for(Long64_t entry = startEntry2; entry < maxEntry2; entry++){    
    cPhotoLdb->GetEntry(entry);
    
    if(TMath::Abs(deltaPhiDegPhotoLdb) < phiCut){

      hDeltaPhiDegPhotoLdb->Fill(deltaPhiDegPhotoLdb);
      hDeltaPhiDegPhotoLdb3->Fill(phiExpected, deltaPhiDegPhotoLdb);

      hDeltaThetaDegPhotoLdb->Fill(deltaThetaDegPhotoLdb);
      hDeltaThetaDegPhotoLdb3->Fill(phiExpected, deltaThetaDegPhotoLdb);
    }
    
    p2.inc(entry, maxEntry2);
  }

  // TObjArray* tArrPhi =  new TObjArray();
  // hDeltaPhiDegPhotoWais2->FitSlicesY(NULL, 1, hDeltaPhiDegPhotoWais2->GetNbinsY(), 0, "QNR", tArrPhi);
  // for(int i=0; i < tArrPhi->GetEntries(); i++){
  //   auto h = (TH1D*)tArrPhi->At(i);
  //   // auto c0 = new TCanvas();
  //   // TString opt = i==0? "" : "same";
  //   // h->Draw(opt);
  //   if(i==2){
  //   }
  // }
  
  // TObjArray* tArrTheta =  new TObjArray();
  // hDeltaThetaDegPhotoWais2->FitSlicesY(NULL, 1, hDeltaThetaDegPhotoWais2->GetNbinsY(), 0, "QNR", tArrTheta);
  // for(int i=0; i < tArrTheta->GetEntries(); i++){
  //   auto h = (TH1D*)tArrTheta->At(i);
  //   // auto c0 = new TCanvas();    
  //   // TString opt = i==0? "" : "same";
  //   // h->Draw(opt);
  //   if(i==2){
  //   }
  // }

  auto c0 = new TCanvas();
  hDeltaPhiDegPhotoLdb3->Draw("colz");
  auto c0a = new TCanvas();
  hDeltaThetaDegPhotoLdb3->Draw("colz");
  
  auto c0b = new TCanvas();
  hDeltaPhiDegPhotoWais3->Draw("colz");
  auto c0c = new TCanvas();
  hDeltaThetaDegPhotoWais3->Draw("colz");
    
  
  // auto c1 = new TCanvas();
  // c1->SetLogy(1);
  hDeltaPhiDegPhotoWais->SetLineColor(kRed);
  hDeltaPhiDegPhotoLdb->SetLineColor(kBlue);
  hDeltaThetaDegPhotoWais->SetLineColor(kMagenta);
  hDeltaThetaDegPhotoLdb->SetLineColor(kCyan);

  hDeltaPhiDegPhotoWais->SetMarkerSize(0);
  hDeltaPhiDegPhotoLdb->SetMarkerSize(0);
  hDeltaThetaDegPhotoWais->SetMarkerSize(0);
  hDeltaThetaDegPhotoLdb->SetMarkerSize(0);


  hDeltaPhiDegPhotoWais->Scale(1./hDeltaPhiDegPhotoWais->Integral());
  hDeltaPhiDegPhotoLdb->Scale(1./hDeltaPhiDegPhotoLdb->Integral());
  hDeltaThetaDegPhotoWais->Scale(1./hDeltaThetaDegPhotoWais->Integral());
  hDeltaThetaDegPhotoLdb->Scale(1./hDeltaThetaDegPhotoLdb->Integral());

  const int nh = 4;
  TH1D* hs[nh] = {hDeltaPhiDegPhotoLdb, hDeltaThetaDegPhotoLdb, hDeltaPhiDegPhotoWais, hDeltaThetaDegPhotoWais};
  auto c1 = RootTools::drawHistsWithStatsBoxes(nh, hs, "", "mre");
  auto l1 = c1->BuildLegend();
  l1->SetNColumns(2);
  l1->Draw();
  hs[0]->SetMaximum(0.13);
  
  hDeltaPhiDegPhotoLdb->SetTitle("ANITA Angular Resolution Using Photogrammetry; #delta#theta or #delta#phi (Degrees); Fraction of events per bin");

}
