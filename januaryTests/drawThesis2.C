#include "TChain.h"
#include "RootTools.h"

void drawThesis2(){

  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_11-57-40.root");
  // TFile* f = TFile::Open("testing/comboStudies/generateAngularResolutionTreePlots_352_2016-03-02_13-05-06.root");
  // TTree* t = (TTree*) f->Get("angResTree");

  {
    TChain* t = new TChain("angResTree");
    t->Add("testing/comboStudies/generateAngularResolutionTreePlots_*");  

    const int numBins = 64; //128;
    const double maxDeg = 2;
    TH1D* hTheta = new TH1D("hTheta", "#delta#theta Distribution of WAIS pulses; #delta#theta (Degrees); Events per bin", numBins, -maxDeg, maxDeg);

    TH1D* hPhi = new TH1D("hPhi", "#delta#phi Distribution of WAIS pulses; #delta#phi (Degrees); Events per bin", numBins, -maxDeg, maxDeg); 

    new TCanvas();
    t->Draw("deltaThetaDeg3:phiExpected", "TMath::Abs(deltaPhiDeg3) < 3 && TMath::Abs(deltaThetaDeg3) < 5", "colz");
    new TCanvas();    
    t->Draw("deltaPhiDeg5:phiExpected", "TMath::Abs(deltaPhiDeg5) < 5", "colz");
    
    // t->Draw("deltaPhiDeg>>h2", "TMath::Abs(deltaPhiDeg) < 5", "colz");
    t->Draw("deltaThetaDeg5>>hTheta", "TMath::Abs(deltaPhiDeg5) < 5", "goff");
    t->Draw("deltaPhiDeg5>>hPhi", "TMath::Abs(deltaPhiDeg5) < 5", "goff");  

    TH1D* hs[2] = {hTheta, hPhi};
    hPhi->SetLineColor(kRed);
    hTheta->SetLineColor(kBlue);  
    auto c1 = RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");

    TLegend* l1 = new TLegend(0.8, 0.8, 1, 1);
    l1->AddEntry(hPhi, "#delta#phi (Degrees)", "l");
    l1->AddEntry(hTheta, "#delta#phi (Degrees)", "l");
    hs[0]->SetTitle("Angular resolution of WAIS pulses");
  }

  {
    TChain* t = new TChain("angResTree");
    t->Add("photogrammetryNumbers/generateAngularResolutionTreePlots_*");  
    std::cout << t->GetEntries() << std::endl;
    const int numBins = 64; //128;
    const double maxDeg = 2;
    TH1D* hTheta = new TH1D("hTheta2", "#delta#theta Distribution of WAIS pulses; #delta#theta (Degrees); Events per bin", numBins, -maxDeg, maxDeg);

    TH1D* hPhi = new TH1D("hPhi2", "#delta#phi Distribution of WAIS pulses; #delta#phi (Degrees); Events per bin", numBins, -maxDeg, maxDeg); 

    // t->Draw("deltaPhiDeg>>h2", "TMath::Abs(deltaPhiDeg) < 5", "colz");
    t->Show(0);
    auto c0 = new TCanvas();
    TH2D* hDeltaThetaDeg2 = new TH2D("hDeltaThetaDeg2", "Difference in expected and measured Elevation vs. Payload Azimuth; #phi_{expected} (Degrees); #Delta#theta_{WAIS} (Degrees) ; Events per bin",
				     256, 0, 360,
				     256, -5, 5);
    t->Draw("deltaThetaDeg:phiExpected>>hDeltaThetaDeg2", "TMath::Abs(globalPhiDeg-zoomPhiDeg) < 5 && TMath::Abs(deltaPhiDeg) < 5 && TMath::Abs(deltaThetaDeg) < 5", "colz");
    auto c2 = new TCanvas();
    TH2D* hDeltaPhiDeg2 = new TH2D("hDeltaPhiDeg2", "Difference in expected and measured Elevation vs. Payload Azimuth; #phi_{expected} (Degrees); #Delta#phi_{WAIS} (Degrees) ; Events per bin",
				   3600, 0, 360,
				   30, -5, 5);
				   // 256, 0, 360,
				   // 256, -5, 5);
    
    // t->Draw("deltaPhiDeg:phiExpected>>hDeltaPhiDeg2", "TMath::Abs(triggeredPhiDeg-zoomPhiDeg) < 999 && TMath::Abs(deltaPhiDeg) < 5", "colz");
    t->Draw("triggeredPhiDeg-zoomPhiDeg:triggeredPhiDeg>>hDeltaPhiDeg2", "TMath::Abs(triggeredPhiDeg-zoomPhiDeg) < 999 && TMath::Abs(deltaPhiDeg) < 5", "colz");    
    // t->Draw("deltaPhiDeg:triggeredPhiDeg>>hDeltaPhiDeg2", "TMath::Abs(triggeredPhiDeg-zoomPhiDeg) < 999 && TMath::Abs(deltaPhiDeg) < 5", "colz");
    // t->Draw("triggeredPhiDeg", "TMath::Abs(deltaPhiDeg) < 5", "colz");        

    new TCanvas();
    t->Draw("triggeredPhiDeg - zoomPhiDeg:triggeredPhiDeg", "TMath::Abs(triggeredPhiDeg-zoomPhiDeg) < 1.24 && TMath::Abs(deltaPhiDeg) < 5", "colz");
    
    t->Draw("deltaThetaDeg>>hTheta2", "TMath::Abs(deltaPhiDeg) < 5", "goff");    
    t->Draw("deltaPhiDeg>>hPhi2", "TMath::Abs(deltaPhiDeg) < 5", "goff");  

    TH1D* hs[2] = {hTheta, hPhi};
    hPhi->SetLineColor(kRed);
    hTheta->SetLineColor(kBlue);  
    auto c1 = RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");

    TLegend* l1 = new TLegend(0.8, 0.8, 1, 1);
    l1->AddEntry(hPhi, "#delta#phi (Degrees)", "l");
    l1->AddEntry(hTheta, "#delta#phi (Degrees)", "l");
    hs[0]->SetTitle("Angular resolution of WAIS pulses");
  }
  
}
