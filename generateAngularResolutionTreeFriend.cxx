// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Add a new tree with pitch/roll corrected expected theta/phi distributions
*************************************************************************************************************** */

#include <TFile.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <THnSparse.h>

#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

//const Double_t bestGuessConstPitch = 

const Double_t provisionalBestPitch = 0.3; //0.471588; //-0.3; //-0.58;
const Double_t provisionalBestRoll = 0.01; //-0.0564042; //-0.01; //0.86;
const Double_t provisionalBestHeading = -0.31; //-0.353051; //0.31; //0.02;

// 0.471588	-0.0564042	-0.353051

int main(int argc, char *argv[])
{

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;    
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* headChain = new TChain("headTree");  
  TChain* angResChain = new TChain("angResTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    angResChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    
  }

  auto fHist = TFile::Open("fitPitchRollOffsetsPlots.root");
  auto hDeltaThetaDeg_pfx = (TProfile*) fHist->Get("hDeltaThetaDeg_pfx");
  auto hDeltaPhiDeg_pfx = (TProfile*) fHist->Get("hDeltaPhiDeg_pfx");
    

  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  
  UInt_t eventNumber = 0;
  angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);    
  angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);
  angResChain->SetBranchAddress("thetaExpected", &thetaExpected);    
  angResChain->SetBranchAddress("phiExpected", &phiExpected);
  
  angResChain->SetBranchAddress("eventNumber", &eventNumber);    

  
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  headChain->BuildIndex("header->eventNumber");  

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");  
  
  TString outFileName = TString::Format("%s_run%d-%dPlots.root", argv[0], firstRun, lastRun);
  TFile* outFile = new TFile(outFileName, "recreate");
  TTree* angResTree = new TTree("angResTreeWithDynamicPitchRoll", "angResTreeWithDynamicPitchRoll");

  Double_t thetaExpectedWithDynamicPitchRoll = 0;
  Double_t phiExpectedWithDynamicPitchRoll = 0;
  Double_t deltaThetaDegWithDynamicPitchRoll = 0;
  Double_t deltaPhiDegWithDynamicPitchRoll = 0;

  Double_t thetaExpectedWithProvisionalBestConstantOffsets = 0;
  Double_t phiExpectedWithProvisionalBestConstantOffsets = 0;
  Double_t deltaThetaDegWithProvisionalBestConstantOffsets = 0;
  Double_t deltaPhiDegWithProvisionalBestConstantOffsets = 0;

  Double_t thetaExpectedBestPossible = 0;
  Double_t phiExpectedBestPossible = 0;
  Double_t deltaThetaDegBestPossible = 0;
  Double_t deltaPhiDegBestPossible = 0;
  
  
  // UInt_t eventNumber = 0;
  Double_t pitch = 0;
  Double_t roll = 0;
  // std::vector<Double_t>* deltaPhiDeg = NULL;

  angResTree->Branch("deltaPhiDegWithDynamicPitchRoll", &deltaPhiDegWithDynamicPitchRoll);
  angResTree->Branch("deltaThetaDegWithDynamicPitchRoll", &deltaThetaDegWithDynamicPitchRoll);  
  angResTree->Branch("thetaExpectedWithDynamicPitchRoll", &thetaExpectedWithDynamicPitchRoll);
  angResTree->Branch("phiExpectedWithDynamicPitchRoll", &phiExpectedWithDynamicPitchRoll);

  angResTree->Branch("deltaPhiDegWithProvisionalBestConstantOffsets", &deltaPhiDegWithProvisionalBestConstantOffsets);
  angResTree->Branch("deltaThetaDegWithProvisionalBestConstantOffsets", &deltaThetaDegWithProvisionalBestConstantOffsets);  
  angResTree->Branch("thetaExpectedWithProvisionalBestConstantOffsets", &thetaExpectedWithProvisionalBestConstantOffsets);
  angResTree->Branch("phiExpectedWithProvisionalBestConstantOffsets", &phiExpectedWithProvisionalBestConstantOffsets);

  angResTree->Branch("deltaPhiDegBestPossible", &deltaPhiDegBestPossible);
  angResTree->Branch("deltaThetaDegBestPossible", &deltaThetaDegBestPossible);  
  angResTree->Branch("thetaExpectedBestPossible", &thetaExpectedBestPossible);
  angResTree->Branch("phiExpectedBestPossible", &phiExpectedBestPossible);

  
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("pitch", &pitch);
  angResTree->Branch("roll", &roll);

  Long64_t nEntries = angResChain->GetEntries();
  Long64_t maxEntry = 0; //5000; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    angResChain->GetEntry(entry);
    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    UsefulAdu5Pat usefulPat(pat);
    pitch = pat->pitch;
    roll = pat->roll;

    usefulPat.pitch += (pitch-0.5);
    usefulPat.roll += (roll+0.57);
    // usefulPat.pitch -= (pitch-0.5);
    // usefulPat.roll += (roll+0.57);

    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpectedWithDynamicPitchRoll,
					   phiExpectedWithDynamicPitchRoll);
    phiExpectedWithDynamicPitchRoll*=TMath::RadToDeg();
    thetaExpectedWithDynamicPitchRoll*=-1*TMath::RadToDeg();

    deltaPhiDegWithDynamicPitchRoll = RootTools::getDeltaAngleDeg(phiExpectedWithDynamicPitchRoll, zoomPhiDeg);
    deltaThetaDegWithDynamicPitchRoll = RootTools::getDeltaAngleDeg(thetaExpectedWithDynamicPitchRoll, zoomThetaDeg);

    // usefulPat.pitch = -provisionalBestPitch;
    // usefulPat.roll = -provisionalBestRoll;
    // usefulPat.heading -= provisionalBestHeading;
    usefulPat.pitch = provisionalBestPitch;
    usefulPat.roll = provisionalBestRoll;
    usefulPat.heading += provisionalBestHeading;

    if(usefulPat.heading < 0) usefulPat.heading += 360;
    if(usefulPat.heading >= 360) usefulPat.heading -= 360;    

    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpectedWithProvisionalBestConstantOffsets,
					   phiExpectedWithProvisionalBestConstantOffsets);
    phiExpectedWithProvisionalBestConstantOffsets*=TMath::RadToDeg();
    thetaExpectedWithProvisionalBestConstantOffsets*=-1*TMath::RadToDeg();

    deltaPhiDegWithProvisionalBestConstantOffsets = RootTools::getDeltaAngleDeg(phiExpectedWithProvisionalBestConstantOffsets, zoomPhiDeg);
    deltaThetaDegWithProvisionalBestConstantOffsets = RootTools::getDeltaAngleDeg(thetaExpectedWithProvisionalBestConstantOffsets, zoomThetaDeg);

    Int_t bin = hDeltaThetaDeg_pfx->GetXaxis()->FindBin(zoomPhiDeg);
    Double_t deltaPhi = hDeltaPhiDeg_pfx->GetBinContent(bin);
    deltaPhiDegBestPossible = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg);
    deltaPhiDegBestPossible -= deltaPhi;
    // if(deltaPhiDegBestPossible < 0) deltaPhiDegBestPossible += 360;
    // if(deltaPhiDegBestPossible >= 360) deltaPhiDegBestPossible -= 360;

    Double_t deltaTheta = hDeltaThetaDeg_pfx->GetBinContent(bin);
    deltaThetaDegBestPossible = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg);
    deltaThetaDegBestPossible -= deltaTheta;
    // if(deltaThetaDegBestPossible < 0) deltaThetaDegBestPossible += 360;
    // if(deltaThetaDegBestPossible >= 360) deltaThetaDegBestPossible -= 360;
    
    
    angResTree->Fill();
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
