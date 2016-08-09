#include "RootTools.h"
#include "ProgressBar.h"
#include "TProfile.h"

void plotZoomPhiDiff(){

  TChain* cPhoto = new TChain("angResTree");
  TChain* cFiltOn = new TChain("angResTree");
  TChain* cFiltOff = new TChain("angResTree");

  cPhoto->Add("filterOffs/photogrammetry/generateAngularResolutionTreePlots_*");  
  cFiltOn->Add("filterOns/phaseCenter/generateAngularResolutionTreePlots_*");
  cFiltOff->Add("filterOffs/phaseCenter/generateAngularResolutionTreePlots_*");
  TString polText = "HPol";
  TString source = "WAIS";

  // cPhoto->Add("filterOffs/photogrammetry/generateAngularResolutionTreeVPOLPlots_*");  
  // cFiltOn->Add("extremeFilter/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
  // cFiltOff->Add("filterOffs/phaseCenter/generateAngularResolutionTreeVPOLPlots_*");
  // TString polText = "VPol";
  // TString source = "LDB";
  
  Double_t zoomPhiDegOn = 0;
  cFiltOn->SetBranchAddress("zoomPhiDeg", &zoomPhiDegOn);
  Double_t zoomPhiDegOff = 0;
  cFiltOff->SetBranchAddress("zoomPhiDeg", &zoomPhiDegOff);

  Double_t deltaPhiDegOn = 0;
  cFiltOn->SetBranchAddress("deltaPhiDeg", &deltaPhiDegOn);
  Double_t deltaPhiDegOff = 0;
  cFiltOff->SetBranchAddress("deltaPhiDeg", &deltaPhiDegOff);  
  Double_t deltaPhiDegPhoto = 0;
  cPhoto->SetBranchAddress("deltaPhiDeg", &deltaPhiDegPhoto);

  Double_t zoomThetaDegOn = 0;
  cFiltOn->SetBranchAddress("zoomThetaDeg", &zoomThetaDegOn);
  Double_t zoomThetaDegOff = 0;
  cFiltOff->SetBranchAddress("zoomThetaDeg", &zoomThetaDegOff);


  Double_t deltaThetaDegOn = 0;
  cFiltOn->SetBranchAddress("deltaThetaDeg", &deltaThetaDegOn);
  Double_t deltaThetaDegOff = 0;
  cFiltOff->SetBranchAddress("deltaThetaDeg", &deltaThetaDegOff);  
  Double_t deltaThetaDegPhoto = 0;
  cPhoto->SetBranchAddress("deltaThetaDeg", &deltaThetaDegPhoto);  

  Double_t phiExpected = 0;
  cFiltOn->SetBranchAddress("phiExpected", &phiExpected);
  Double_t thetaExpected = 0;
  cFiltOn->SetBranchAddress("thetaExpected", &thetaExpected);
  
  Double_t heading = 0;
  cFiltOn->SetBranchAddress("heading", &heading);
  Double_t realTime = 0;
  cFiltOn->SetBranchAddress("realTime", &realTime);

  Long64_t nEntries = cFiltOn->GetEntries();
  Long64_t maxEntry = 0; //5; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  double dPhi = 0.05;  
  double phiRange = 6;

  double dTheta = 0.05;  
  double thetaRange = 6;

  const int numBinsAzimuth = 256;
  
  TH1D* hFilterPhiDiff = new TH1D("hfilterPhiDiff", "#phi_{Filtered} - #phi_{Unfiltered}; #delta#phi (Degrees); Events per bin", phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH2D* hFilterPhiDiff2 = new TH2D("hfilterPhiDiff2",
				   polText + " #delta#phi_{"+source+"} (Filtered); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				   numBinsAzimuth, 0, 360,
				   phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH2D* hFilterPhiDiff2a = new TH2D("hfilterPhiDiff2a",
				   polText + " #delta#phi_{"+source+"} (Unfiltered); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				   numBinsAzimuth, 0, 360,
				   phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);  


  TH2D* hFilterPhiDiff2b = new TH2D("hfilterPhiDiff2b",
				    polText + " #delta#phi_{"+source+"} (Photogrammetry); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				    numBinsAzimuth, 0, 360,
				    phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);  
 
  TH2D* hFilterPhiDiff3 = new TH2D("hfilterPhiDiff3",
				   "#phi_{Filtered} - #phi_{Unfiltered}; Heading (Degrees); #delta#phi (Degrees); Events per bin",
				   numBinsAzimuth, 0, 360,
				   phiRange/dPhi, -phiRange/2 - dPhi/2, phiRange/2 - dPhi/2);

  TH1D* hFilterThetaDiff = new TH1D("hfilterThetaDiff", "#theta_{Filtered} - #theta_{Unfiltered}; #delta#theta (Degrees); Events per bin", thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  TH2D* hFilterThetaDiff2 = new TH2D("hfilterThetaDiff2",
				   polText + " #delta#phi_{"+source+"} (Filtered); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				   numBinsAzimuth, 0, 360,
				   thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  TH2D* hFilterThetaDiff2a = new TH2D("hfilterThetaDiff2a",
				   polText + " #delta#phi_{"+source+"} (Unfiltered); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				      numBinsAzimuth, 0, 360,
				      thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);


  TH2D* hFilterThetaDiff2b = new TH2D("hfilterThetaDiff2b",
				      polText + " #delta#phi_{"+source+"} (Photogrammery); #phi_{expected} (Degrees); #delta#phi (Degrees); Events per bin",
				      numBinsAzimuth, 0, 360,
				      thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  TH2D* hFilterThetaDiff3 = new TH2D("hfilterThetaDiff3",
				   "#theta_{Filtered} - #theta_{Unfiltered}; Heading (Degrees); #delta#theta (Degrees); Events per bin",
				   numBinsAzimuth, 0, 360,
				   thetaRange/dTheta, -thetaRange/2 - dTheta/2, thetaRange/2 - dTheta/2);

  const double maxDeltaPhiDeg = 3;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    cFiltOn->GetEntry(entry);
    cFiltOff->GetEntry(entry);
    cPhoto->GetEntry(entry);
    Double_t filterPhiDiff = RootTools::getDeltaAngleDeg(zoomPhiDegOn, zoomPhiDegOff);
    Double_t filterThetaDiff = RootTools::getDeltaAngleDeg(zoomThetaDegOn, zoomThetaDegOff);

    if(deltaPhiDegOff < maxDeltaPhiDeg){
    
      hFilterPhiDiff->Fill(filterPhiDiff);
      // hFilterPhiDiff2->Fill(thetaExpected, deltaPhiDegOn);
      // hFilterPhiDiff2a->Fill(thetaExpected, deltaPhiDegOff);
      // hFilterPhiDiff2b->Fill(thetaExpected, deltaPhiDegPhoto);
      hFilterPhiDiff2->Fill(phiExpected, deltaPhiDegOn);
      hFilterPhiDiff2a->Fill(phiExpected, deltaPhiDegOff);
      hFilterPhiDiff2b->Fill(phiExpected, deltaPhiDegPhoto);
      

      hFilterPhiDiff3->Fill(heading, filterPhiDiff);
      // hFilterPhiDiff2->Fill(phiExpected, filterPhiDiff);
      // hFilterPhiDiff3->Fill(heading, filterPhiDiff);

      hFilterThetaDiff->Fill(filterThetaDiff);
      // hFilterThetaDiff2->Fill(thetaExpected, deltaThetaDegOn);
      // hFilterThetaDiff2a->Fill(thetaExpected, deltaThetaDegOff);
      // hFilterThetaDiff2b->Fill(thetaExpected, deltaThetaDegPhoto);
      hFilterThetaDiff2->Fill(phiExpected, deltaThetaDegOn);
      hFilterThetaDiff2a->Fill(phiExpected, deltaThetaDegOff);
      hFilterThetaDiff2b->Fill(phiExpected, deltaThetaDegPhoto);      
      
      hFilterThetaDiff3->Fill(heading, filterThetaDiff);
      // hFilterThetaDiff2->Fill(phiExpected, filterThetaDiff);
      // hFilterThetaDiff3->Fill(heading, filterThetaDiff);
    }
    
    p.inc(entry, maxEntry);
  }

  // c1->SetLogy(1);
  hFilterPhiDiff->SetLineColor(kRed);
  hFilterThetaDiff->SetLineColor(kBlue);
  hFilterPhiDiff->SetMarkerSize(0);
  hFilterThetaDiff->SetMarkerSize(0);

  const int nProjs = 2;
  TH1D* hProjs[nProjs] = {hFilterThetaDiff, hFilterPhiDiff};
  auto c1 = RootTools::drawHistsWithStatsBoxes(nProjs, hProjs, "", "mre");
  auto l1 = c1->BuildLegend();
  l1->Draw();
  hFilterThetaDiff->SetTitle(polText + " Change in Interferometric Peak when Filtering; Angular Shift (Degrees); Events per bin");

  auto c2 = new TCanvas();
  // c2->SetLogz(1);
  hFilterPhiDiff2->Draw("colz");

  auto c2a = new TCanvas();
  // c2->SetLogz(1);
  hFilterPhiDiff2a->Draw("colz");


  TProfile* x = hFilterPhiDiff2->ProfileX();
  TProfile* x2 = hFilterPhiDiff2a->ProfileX();
  TProfile* x3 = hFilterPhiDiff2b->ProfileX();
  x->SetTitle(polText + " #delta#phi_{LDB} Filtered");
  x2->SetTitle(polText + " #delta#phi_{LDB} Unfiltered");
  x3->SetTitle(polText + " #delta#phi_{LDB} Photogrammetry");  

  x->Rebin(2);
  x2->Rebin(2);
  x3->Rebin(2);

  // x->Rebin(8);
  // x2->Rebin(8);


  // c1->SetLogy(1);
  x->SetLineColor(kRed);
  x2->SetLineColor(kBlue);
  x3->SetLineColor(kMagenta);  
  x->SetMarkerSize(0);
  x2->SetMarkerSize(0);
  x3->SetMarkerSize(0);  
  
  const int nF = 3;
  TProfile* hFs[nF] = {x2, x, x3};
  x2->SetMaximum(1.8);
  x2->SetMinimum(-1.8);
  auto c2b = RootTools::drawHistsWithStatsBoxes(nF, (TH1D**)hFs, "", "mre");
  auto l2b = c2b->BuildLegend();
  l2b->Draw();
  // x->SetTitle(polText + " Change in Interferometric Peak when Filtering; Angular Shift (Degrees); Events per bin");
  x2->SetTitle(polText + " Comparison of #delta#phi_{"+source+"}; #phi_{LDB};  #phi_{"+source+"}(Degree);");  

  auto c3 = new TCanvas();
  // c3->SetLogz(1);
  hFilterPhiDiff3->Draw("colz");

  auto c3a = new TCanvas();
  // c2->SetLogz(1);
  hFilterThetaDiff2a->Draw("colz");


  TProfile* y = hFilterThetaDiff2->ProfileX();
  TProfile* y2 = hFilterThetaDiff2a->ProfileX();
  TProfile* y3 = hFilterThetaDiff2b->ProfileX();
  y->SetTitle(polText + " #delta#theta_{LDB} Filtered");
  y2->SetTitle(polText + " #delta#theta_{LDB} Unfiltered");
  y3->SetTitle(polText + " #delta#theta_{LDB} Photogrammetry");    
  y->Rebin(2);
  y2->Rebin(2);
  y3->Rebin(2);
    
  // x->Rebin(8);
  // x2->Rebin(8);


  // c1->SetLogy(1);
  y->SetLineColor(kRed);
  y2->SetLineColor(kBlue);
  y3->SetLineColor(kMagenta);  
  y->SetMarkerSize(0);
  y2->SetMarkerSize(0);
  y3->SetMarkerSize(0);  
  
  const int nF2 = 3;
  TProfile* hF2s[nF2] = {y2, y, y3};

  if(polText.Contains("H")){
    y2->SetMaximum(0.6);
    y2->SetMinimum(-0.6);
  }
  else{
    y2->SetMaximum(0.2);
    y2->SetMinimum(-1.8);
  }
  
  auto c3b = RootTools::drawHistsWithStatsBoxes(nF2, (TH1D**)hF2s, "", "mre");
  auto l3b = c3b->BuildLegend();
  l3b->Draw();
  // x->SetTitle(polText + " Change in Interferometric Peak when Filtering; Angular Shift (Degrees); Events per bin");
  y2->SetTitle(polText + " Comparison of #delta#theta_{"+source+"}; #theta_{LDB};  #theta_{"+source+"}(Degree);");

  // auto c4 = new TCanvas();
  // // c4->SetLogy(1);

  auto c5 = new TCanvas();
  // c5->SetLogz(1);
  hFilterThetaDiff2->Draw("colz");  

  auto c6 = new TCanvas();
  // c6->SetLogz(1);
  hFilterThetaDiff3->Draw("colz");  
  
}
