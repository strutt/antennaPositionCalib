// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Loop over the individual events with a function in MINUIT. Note that this program assumes no constant pitch roll offsets in ANITA framework. 
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
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>


#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

//const Double_t bestGuessConstPitch = 

const Double_t provisionalBestPitch = -0.3; //-0.58;
const Double_t provisionalBestRoll = -0.01; //0.86;
const Double_t provisionalBestHeading = 0.31; //0.02;


Double_t funcToMin(const Double_t *vars);
Long64_t maxEntry = 0;
Long64_t startEntry = 0;

Double_t zoomPhiDeg = 0;
Double_t zoomThetaDeg = 0;
UInt_t eventNumber = 0;
RawAnitaHeader* header = NULL;
Adu5Pat* pat = NULL;


TChain* gpsChain = NULL;
TChain* headChain = NULL;
TChain* angResChain = NULL;

const Double_t deltaPhiMax = 2;
const Double_t deltaThetaMax = 2;

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

  gpsChain = new TChain("adu5PatTree");
  headChain = new TChain("headTree");  
  angResChain = new TChain("angResTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run);
    angResChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    
  }    

  zoomPhiDeg = 0;
  zoomThetaDeg = 0;
  eventNumber = 0;
  angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
  angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);
  angResChain->SetBranchAddress("eventNumber", &eventNumber);

  
  header = NULL;
  headChain->SetBranchAddress("header", &header);
  headChain->BuildIndex("header->eventNumber");  

  pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");
    
  Long64_t nEntries = angResChain->GetEntries();
  maxEntry = 0; //5000; //5000;
  startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;



  // Mostly copied from the ROOT Minuit tutorial 

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  
  Int_t numVars = 3;
  ROOT::Math::Functor FuncToMin(&funcToMin, numVars);

  Double_t stepSize = 1e-3;
  std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);
  
  // starting point
  std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);
  variables.at(0) = provisionalBestPitch;
  variables.at(1) = provisionalBestRoll;
  variables.at(2) = provisionalBestHeading;

  
  min->SetFunction(FuncToMin);
  
  min->SetVariable(0, "pitch offset", variables[0], step[0]);
  min->SetVariable(1, "roll offset", variables[1], step[1]);  
  min->SetVariable(2, "heading offset", variables[2], step[2]);  

  // do the minimization
  min->Minimize(); 

  
  return 0;
  
}


Double_t funcToMin(const Double_t *vars){

  Double_t pitchOffset = vars[0];
  Double_t rollOffset = vars[1];
  Double_t headingOffset = vars[2];

  Double_t sumDiffSq = 0;

  // ProgressBar p(maxEntry-startEntry);
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){    
    angResChain->GetEntry(entry);
    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    UsefulAdu5Pat usefulPat(pat);

    usefulPat.heading += headingOffset;
    if(usefulPat.heading >= 360) usefulPat.heading -= 360;
    if(usefulPat.heading < 0) usefulPat.heading += 360;    
    usefulPat.pitch += pitchOffset;
    usefulPat.roll += rollOffset;
    
    Double_t thetaExpected = 0;
    Double_t phiExpected = 0;    
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
    thetaExpected *= -1*TMath::RadToDeg();
    phiExpected *= TMath::RadToDeg();    
    
    Double_t deltaPhi = RootTools::getDeltaAngleDeg(zoomPhiDeg, phiExpected);
    Double_t deltaTheta = RootTools::getDeltaAngleDeg(zoomThetaDeg, thetaExpected);

    // if(TMath::Abs(deltaPhi) > 3 || TMath::Abs(deltaTheta) > 3){
    if(TMath::Abs(deltaPhi) > deltaPhiMax || TMath::Abs(deltaTheta) > deltaThetaMax){
      continue;
    }
    
    Double_t diffSq = deltaPhi*deltaPhi + deltaTheta*deltaTheta;
    // std::cout << deltaPhi << "\t" << deltaTheta << "\t" << diffSq << std::endl;
    
    sumDiffSq += diffSq;

    // std::cout << sumDiffSq/(maxEntry-startEntry) << std::endl;
    
    // p++;
  }

  std::cout << pitchOffset << "\t" << rollOffset << "\t" << headingOffset << "\t"
	    << sumDiffSq << std::endl;
  return sumDiffSq;

}
