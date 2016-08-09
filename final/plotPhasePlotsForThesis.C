#include "RootTools.h"
#include "ProgressBar.h"

void plotPhasePlotsForThesis(){

  TChain* cPhaseLdb = new TChain("angResTree");
  TChain* cPhaseWais = new TChain("angResTree");

  cPhaseWais->Add("filterOffs/phaseCenter/generateAngularResolutionTreePlots_*");
  cPhaseLdb->Add("filterOffs/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
    
  Double_t deltaPhiDegPhaseWais = 0;  
  cPhaseWais->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhaseWais);
  Double_t deltaPhiDegPhaseLdb = 0;  
  cPhaseLdb->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhaseLdb);

  Double_t deltaThetaDegPhaseWais = 0;  
  cPhaseWais->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhaseWais);
  Double_t deltaThetaDegPhaseLdb = 0;  
  cPhaseLdb->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhaseLdb);  
  
  
  
  Double_t phiExpected = 0;
  cPhaseWais->SetBranchAddress("phiExpected", &phiExpected);
  Double_t heading = 0;
  cPhaseWais->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  cPhaseWais->SetBranchAddress("realTime", &realTime);

  Long64_t nEntries = cPhaseWais->GetEntries();
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

  
  TH1D* hDeltaPhiDegPhaseWais = new TH1D("hDeltaPhiDegPhaseWais", "HPol #delta#phi_{WAIS}", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);
  TH2D* hDeltaPhiDegPhaseWais3 = new TH2D("hDeltaPhiDegPhaseWais3", "HPol #phi_{WAIS} - #phi_{measured} Fitted Phase Centre Antenna Positions; #delta#phi_{WAIS} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH1D* hDeltaThetaDegPhaseWais = new TH1D("hDeltaThetaDegPhaseWais", "HPol #delta#theta_{WAIS}", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhaseWais3 = new TH2D("hDeltaThetaDegPhaseWais3", "HPol #theta_{WAIS} - #theta_{measured} Fitted Phase Centre Antenna Positions; #delta#theta_{WAIS} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);


  
  TH1D* hDeltaPhiDegPhaseLdb = new TH1D("hDeltaPhiDegPhaseLdb", "VPol #delta#phi_{LDB}", phiRange/(dPhi*rebin), -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);  
  TH2D* hDeltaPhiDegPhaseLdb3 = new TH2D("hDeltaPhiDegPhaseLdb3", "VPol #phi_{LDB} - #phi_{measured} Fitted Phase Centre Antenna Positions; #delta#phi_{LDB} (Degrees); #delta#phi (Degrees); Events per bin", nBinsPhi, 0, 360, phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);


  TH1D* hDeltaThetaDegPhaseLdb = new TH1D("hDeltaThetaDegPhaseLdb", "VPol #delta#theta_{LDB}", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  TH2D* hDeltaThetaDegPhaseLdb3 = new TH2D("hDeltaThetaDegPhaseLdb3", "VPol #theta_{LDB} - #theta_{measured} Fitted Phase CentreAntenna Positions; #delta#theta_{LDB} (Degrees); #delta#theta (Degrees); Events per bin", nBinsPhi, 0, 360, thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);
  
  const double phiCut = 3;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cPhaseWais->GetEntry(entry);

    if(TMath::Abs(deltaPhiDegPhaseWais) < phiCut){

      hDeltaPhiDegPhaseWais->Fill(deltaPhiDegPhaseWais);
      hDeltaPhiDegPhaseWais3->Fill(phiExpected, deltaPhiDegPhaseWais);

      hDeltaThetaDegPhaseWais->Fill(deltaThetaDegPhaseWais);      
      hDeltaThetaDegPhaseWais3->Fill(phiExpected, deltaThetaDegPhaseWais);
      
    }
    p.inc(entry, maxEntry);    
  }

  phiExpected = 0;
  cPhaseLdb->SetBranchAddress("phiExpected", &phiExpected);  

  Long64_t nEntries2 = cPhaseLdb->GetEntries();
  Long64_t maxEntry2 = 0; //5; //1000; //10000;
  Long64_t startEntry2 = 0;
  if(maxEntry2<=0 || maxEntry2 > nEntries2) maxEntry2 = nEntries2;
  std::cout << "Processing " << maxEntry2 << " of " << nEntries2 << " entries." << std::endl;
  ProgressBar p2(maxEntry2-startEntry2);
  
  for(Long64_t entry = startEntry2; entry < maxEntry2; entry++){    
    cPhaseLdb->GetEntry(entry);
    
    if(TMath::Abs(deltaPhiDegPhaseLdb) < phiCut){

      hDeltaPhiDegPhaseLdb->Fill(deltaPhiDegPhaseLdb);
      hDeltaPhiDegPhaseLdb3->Fill(phiExpected, deltaPhiDegPhaseLdb);

      hDeltaThetaDegPhaseLdb->Fill(deltaThetaDegPhaseLdb);
      hDeltaThetaDegPhaseLdb3->Fill(phiExpected, deltaThetaDegPhaseLdb);
    }
    
    p2.inc(entry, maxEntry2);
  }

  // TObjArray* tArrPhi =  new TObjArray();
  // hDeltaPhiDegPhaseWais2->FitSlicesY(NULL, 1, hDeltaPhiDegPhaseWais2->GetNbinsY(), 0, "QNR", tArrPhi);
  // for(int i=0; i < tArrPhi->GetEntries(); i++){
  //   auto h = (TH1D*)tArrPhi->At(i);
  //   // auto c0 = new TCanvas();
  //   // TString opt = i==0? "" : "same";
  //   // h->Draw(opt);
  //   if(i==2){
  //   }
  // }
  
  // TObjArray* tArrTheta =  new TObjArray();
  // hDeltaThetaDegPhaseWais2->FitSlicesY(NULL, 1, hDeltaThetaDegPhaseWais2->GetNbinsY(), 0, "QNR", tArrTheta);
  // for(int i=0; i < tArrTheta->GetEntries(); i++){
  //   auto h = (TH1D*)tArrTheta->At(i);
  //   // auto c0 = new TCanvas();    
  //   // TString opt = i==0? "" : "same";
  //   // h->Draw(opt);
  //   if(i==2){
  //   }
  // }

  auto c0 = new TCanvas();
  hDeltaPhiDegPhaseLdb3->Draw("colz");
  auto c0a = new TCanvas();
  hDeltaThetaDegPhaseLdb3->Draw("colz");
  
  auto c0b = new TCanvas();
  hDeltaPhiDegPhaseWais3->Draw("colz");
  auto c0c = new TCanvas();
  hDeltaThetaDegPhaseWais3->Draw("colz");
    
  
  // auto c1 = new TCanvas();
  // c1->SetLogy(1);
  hDeltaPhiDegPhaseWais->SetLineColor(kRed);
  hDeltaPhiDegPhaseLdb->SetLineColor(kBlue);
  hDeltaThetaDegPhaseWais->SetLineColor(kMagenta);
  hDeltaThetaDegPhaseLdb->SetLineColor(kCyan);

  hDeltaPhiDegPhaseWais->SetMarkerSize(0);
  hDeltaPhiDegPhaseLdb->SetMarkerSize(0);
  hDeltaThetaDegPhaseWais->SetMarkerSize(0);
  hDeltaThetaDegPhaseLdb->SetMarkerSize(0);


  hDeltaPhiDegPhaseWais->Scale(1./hDeltaPhiDegPhaseWais->Integral());
  hDeltaPhiDegPhaseLdb->Scale(1./hDeltaPhiDegPhaseLdb->Integral());
  hDeltaThetaDegPhaseWais->Scale(1./hDeltaThetaDegPhaseWais->Integral());
  hDeltaThetaDegPhaseLdb->Scale(1./hDeltaThetaDegPhaseLdb->Integral());

  const int nh = 4;
  TH1D* hs[nh] = {hDeltaPhiDegPhaseLdb, hDeltaThetaDegPhaseLdb, hDeltaPhiDegPhaseWais, hDeltaThetaDegPhaseWais};
  auto c1 = RootTools::drawHistsWithStatsBoxes(nh, hs, "", "mre");
  auto l1 = c1->BuildLegend();
  // l1->SetNColumns(2);
  l1->Draw();
  hs[0]->SetMaximum(0.13);
  
  hDeltaPhiDegPhaseLdb->SetTitle("ANITA Angular Resolution Using Fitted Phase Centres; #delta#theta or #delta#phi (Degrees); Fraction of events per bin");

}
