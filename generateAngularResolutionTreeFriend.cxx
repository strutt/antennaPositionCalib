// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Add a new tree with pitch/roll corrected expected theta/phi distributions
*****************************************************************************************************************/

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

const Double_t provisionalBestHeading = 0.31; //-0.31; //-0.353051; //0.31; //0.02;
const Double_t provisionalBestPitch = -0.29; //0.3; //0.471588; //-0.3; //-0.58;
const Double_t provisionalBestRoll = -0.02; //0.01; //-0.0564042; //-0.01; //0.86;

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

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnitaIIINumbers(1);
  
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* headChain = new TChain("headTree");
  TChain* angResChain = new TChain("angResTree");
  TChain* deltaTChain = new TChain("deltaTTree");

  for(Int_t run=firstRun; run<=lastRun; run++){

    // Only using event number from this
    TString fileName = TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run);
    deltaTChain->Add(fileName);

    fileName = TString::Format("photogrammetryNoPitchRoll/generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    // fileName = TString::Format("generateAngularResolutionTree_run%d-%dBenPlots.root", run, run);
    angResChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    
  }

  // auto fHist = TFile::Open("fitPitchRollOffsetsPlots.root");
  // auto hDeltaThetaDeg_pfx = (TProfile*) fHist->Get("hDeltaThetaDeg_pfx");
  // auto hDeltaPhiDeg_pfx = (TProfile*) fHist->Get("hDeltaPhiDeg_pfx");
    

  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  
  UInt_t eventNumber = 0;
  angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
  angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);
  angResChain->SetBranchAddress("thetaExpected", &thetaExpected);
  angResChain->SetBranchAddress("phiExpected", &phiExpected);
  deltaTChain->SetBranchAddress("eventNumber", &eventNumber);
  
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  headChain->BuildIndex("header->eventNumber");  

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");  
  
  TString outFileName = TString::Format("%s_run%d-%dPlots.root", argv[0], firstRun, lastRun);
  TFile* outFile = new TFile(outFileName, "recreate");
  TTree* angResTree = new TTree("angResTreeFriend", "angResTreeFriend");

  Double_t thetaExpected0 = 0;
  Double_t phiExpected0 = 0;
  Double_t deltaThetaDeg0 = 0;
  Double_t deltaPhiDeg0 = 0;

  Double_t thetaExpectedTilted = 0;
  Double_t phiExpectedTilted = 0;
  Double_t deltaThetaDegTilted = 0;
  Double_t deltaPhiDegTilted = 0;
  
  // UInt_t eventNumber = 0;
  Double_t pitch = 0;
  Double_t roll = 0;
  // std::vector<Double_t>* deltaPhiDeg = NULL;

  angResTree->Branch("deltaPhiDeg0", &deltaPhiDeg0);
  angResTree->Branch("deltaThetaDeg0", &deltaThetaDeg0);  
  angResTree->Branch("thetaExpected0", &thetaExpected0);
  angResTree->Branch("phiExpected0", &phiExpected0);

  angResTree->Branch("deltaPhiDegTilted", &deltaPhiDegTilted);
  angResTree->Branch("deltaThetaDegTilted", &deltaThetaDegTilted);  
  angResTree->Branch("thetaExpectedTilted", &thetaExpectedTilted);
  angResTree->Branch("phiExpectedTilted", &phiExpectedTilted);

  // angResTree->Branch("deltaPhiDegBestPossible", &deltaPhiDegBestPossible);
  // angResTree->Branch("deltaThetaDegBestPossible", &deltaThetaDegBestPossible);  
  // angResTree->Branch("thetaExpectedBestPossible", &thetaExpectedBestPossible);
  // angResTree->Branch("phiExpectedBestPossible", &phiExpectedBestPossible);

  
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("pitch", &pitch);
  angResTree->Branch("roll", &roll);

  // Long64_t nEntries = angResChain->GetEntries();
  Long64_t nEntries = deltaTChain->GetEntries();  
  Long64_t maxEntry = 0; //5000; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  if(!(AnitaStaticAdu5Offsets::pitch==0 && AnitaStaticAdu5Offsets::roll==0 && AnitaStaticAdu5Offsets::heading==0)){
    std::cerr << "Static offsets: " << std::endl;
    std::cerr << AnitaStaticAdu5Offsets::pitch << "\t" << AnitaStaticAdu5Offsets::roll << "\t" << AnitaStaticAdu5Offsets::heading << "\t" << std::endl;
    std::cerr << "I've made a huge mistake... Please edit the AnitaStaticAdu5Offsets namespace (in UsefulAdu5Pat.h) entries to all be zero." << std::endl;
    return 0;
  }

  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry); 
    angResChain->GetEntry(entry);   
    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    if((header->trigType & 1)==1 && header->eventNumber){;// < 60.95e6 && header->eventNumber > 60.85e6){
    
      UsefulAdu5Pat usefulPat(pat);
      pitch = pat->pitch;
      roll = pat->roll;

      // usefulPat.pitch += (pitch-0.5);
      // usefulPat.roll += (roll+0.57);
      // usefulPat.pitch -= (pitch-0.5);
      // usefulPat.roll += (roll+0.57);

      std::cout << usefulPat.pitch << "\t" << usefulPat.roll << std::endl;

      usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0,
					     phiExpected0);
      phiExpected0*=TMath::RadToDeg();
      thetaExpected0*=-1*TMath::RadToDeg();

      deltaPhiDeg0 = RootTools::getDeltaAngleDeg(phiExpected0, zoomPhiDeg);
      deltaThetaDeg0 = RootTools::getDeltaAngleDeg(thetaExpected0, zoomThetaDeg);


      usefulPat.pitch = provisionalBestPitch;
      usefulPat.roll = provisionalBestRoll;
      usefulPat.heading += provisionalBestHeading;
      if(usefulPat.heading < 0) usefulPat.heading += 360;
      if(usefulPat.heading >= 360) usefulPat.heading -= 360;

      usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpectedTilted,
					     phiExpectedTilted);
      phiExpectedTilted*=TMath::RadToDeg();
      thetaExpectedTilted*=-1*TMath::RadToDeg();

      deltaPhiDegTilted = RootTools::getDeltaAngleDeg(phiExpectedTilted, zoomPhiDeg);
      deltaThetaDegTilted = RootTools::getDeltaAngleDeg(thetaExpectedTilted, zoomThetaDeg);

      // Int_t bin = hDeltaThetaDeg_pfx->GetXaxis()->FindBin(zoomPhiDeg);
      // Double_t deltaPhi = hDeltaPhiDeg_pfx->GetBinContent(bin);
      // deltaPhiDegBestPossible = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg);
      // deltaPhiDegBestPossible -= deltaPhi;
      // // if(deltaPhiDegBestPossible < 0) deltaPhiDegBestPossible += 360;
      // // if(deltaPhiDegBestPossible >= 360) deltaPhiDegBestPossible -= 360;

      // Double_t deltaTheta = hDeltaThetaDeg_pfx->GetBinContent(bin);
      // deltaThetaDegBestPossible = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg);
      // deltaThetaDegBestPossible -= deltaTheta;
      // // if(deltaThetaDegBestPossible < 0) deltaThetaDegBestPossible += 360;
      // // if(deltaThetaDegBestPossible >= 360) deltaThetaDegBestPossible -= 360;
    
    
      angResTree->Fill();
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
