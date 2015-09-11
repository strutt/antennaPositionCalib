// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to find peak cross correlation offsets between antenna pairs in pulses from Wais Divide.
*************************************************************************************************************** */

#include <TFile.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TProfile2D.h>

#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

int main(int argc, char *argv[])
{

  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1] << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);

  TChain* deltaTChain = new TChain("deltaTTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateDeltaTTree_run336-336Plots.root");
    deltaTChain->Add(fileName);
  }

  std::vector<std::vector<Double_t> >* correlationDeltaTs = NULL;
  std::vector<std::vector<Double_t> >* correlationValues = NULL;  
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  std::vector<Double_t>* deltaPhiDeg = NULL;
  deltaTChain->SetBranchAddress("correlationDeltaTs", &correlationDeltaTs);
  deltaTChain->SetBranchAddress("correlationValues", &correlationValues);  
  deltaTChain->SetBranchAddress("thetaExpected", &thetaExpected);
  deltaTChain->SetBranchAddress("phiExpected", &phiExpected);
  deltaTChain->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);  
  
  Long64_t nEntries = deltaTChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  
  // Int_t upsampleFactor = 32;
  // CrossCorrelator* cc = new CrossCorrelator(upsampleFactor);
  // Double_t phiDegMin = 0;
  // Double_t phiDegMax = 360;
  Double_t phiDegMin = 0;
  Double_t phiDegMax = 360;
  const Int_t numBinsPhi = 512;
  Double_t thetaDegMin = -7;
  Double_t thetaDegMax = -5;
  const Int_t numBinsTheta = 45;
  Double_t maxDeltaPhiDeg = 22.5*2;
  
  TProfile2D* prof = new TProfile2D("hProf", "asdasd",
				    numBinsPhi, phiDegMin, phiDegMax,
				    numBinsTheta, thetaDegMin, thetaDegMax);
  TH2D* hCorrDts = new TH2D("hCorrDts", "asd", 1024, -10, 10, 1024, 0, 1);
  
  std::cout << "Generating lookup TProfile2Ds" << std::endl;
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry);
    if(TMath::Abs(deltaPhiDeg->at(0)) < maxDeltaPhiDeg){
      if(correlationValues->at(0).at(6) > 0.4){
	if(correlationDeltaTs->at(0).at(6) > 3 && correlationDeltaTs->at(0).at(6) < 5.5){
	  prof->Fill(phiExpected, thetaExpected, correlationDeltaTs->at(0).at(6));
	  // prof->Fill(deltaPhiDeg->at(0), thetaExpected, correlationDeltaTs->at(0).at(6));
	}
	hCorrDts->Fill(correlationDeltaTs->at(0).at(6), correlationValues->at(0).at(6));
      }
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}


// Double_t getFigureOfMerit(UsefulAdu5Pat* usefulPat,
// 			  Int_t ant1,
// 			  Int_t ant2,
// 			  Double_t deltaR,
// 			  Double_t deltaZ,
// 			  Double_t deltaPhi){



  
  
  
  
// }

