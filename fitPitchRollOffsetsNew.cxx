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


Double_t maxDeltaAngleDeg = 22.5;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;  

UInt_t eventNumber = 0;
Double_t correlationDeltaTs[NUM_COMBOS] = {0};
// Double_t correlationValues[NUM_COMBOS] = {0};  
Double_t correlationDeltaTsClose[NUM_COMBOS] = {0};
// Double_t correlationValuesClose[NUM_COMBOS] = {0};  

Long64_t maxEntry = 0; //3000;

TChain* deltaTChain = NULL;
Adu5Pat* pat = NULL;
TChain* gpsChain = NULL;
RawAnitaHeader* header = NULL;
TChain* headChain = NULL;
AnitaGeomTool* geom = NULL;

Double_t sumOverSquaredDiffs(const Double_t* vars);

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

  geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);

  deltaTChain = new TChain("deltaTTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    deltaTChain->Add(TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run));
  }
  
  deltaTChain->SetBranchAddress("eventNumber", &eventNumber);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTs[%d]", NUM_COMBOS), correlationDeltaTs);
  // deltaTChain->SetBranchAddress(TString::Format("correlationValues[%d]", NUM_COMBOS), correlationValuesClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTsClose[%d]", NUM_COMBOS), correlationDeltaTsClose);
  // deltaTChain->SetBranchAddress(TString::Format("correlationValuesClose[%d]", NUM_COMBOS), correlationValues);
  
  pat =0;
  gpsChain = RootTools::getAdu5PatChain(firstRun, lastRun, pat);
  gpsChain->BuildIndex("realTime");

  header = NULL;
  headChain = RootTools::getHeadChain(firstRun, lastRun, header);
  headChain->BuildIndex("eventNumber");
  
  Long64_t nEntries = deltaTChain->GetEntries();
  maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;

  CrossCorrelator* cc = new CrossCorrelator();
  
  for(Int_t combo=0; combo < NUM_COMBOS; combo++){
    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);

    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
  }
  delete cc;


  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  
  Int_t numVars = 3;
  ROOT::Math::Functor FuncToMin(&sumOverSquaredDiffs, numVars);

  Double_t stepSize = 1e-3;
  std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);
  
  // starting point
  std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);
  variables.at(0) = 0;
  variables.at(1) = 0;
  variables.at(2) = 0;

  
  min->SetFunction(FuncToMin);
  
  min->SetVariable(0, "pitch offset", variables[0], step[0]);
  min->SetVariable(1, "roll offset", variables[1], step[1]);  
  min->SetVariable(2, "heading offset", variables[2], step[2]);  


  // Time it
  TStopwatch watch;
  watch.Start(kTRUE);

  // do the minimization
  min->Minimize(); 

  // Time!
  watch.Start(kFALSE);
  Int_t seconds = Int_t(watch.RealTime());
  Int_t hours = seconds / 3600;
  hours = hours < 0 ? 0 : hours;
  seconds = seconds - hours * 3600;
  Int_t mins = seconds / 60;
  mins = mins < 0 ? 0 : mins;
  seconds = seconds - mins * 60;
  fprintf(stderr, "Minimization took %02d:%02d:%02d\n", hours, mins, seconds);
  
  std::cout << "Minimum = " << min->MinValue() << std::endl;

  return 0;

}

Double_t sumOverSquaredDiffs(const Double_t* vars){

  Double_t pitchOffset = vars[0];
  Double_t rollOffset = vars[1];
  Double_t headingOffset = vars[2];

  Double_t sumOfSquaredDiffs = 0;
  Double_t count = 0;

  for(Long64_t entry = 0; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry);

    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    for(int combo=0; combo<NUM_COMBOS; combo++){
      Int_t ant1=ant1s.at(combo);
      Int_t ant2=ant2s.at(combo);

      if(correlationDeltaTsClose[combo]==correlationDeltaTs[combo]){
	UsefulAdu5Pat usefulPat(pat);
	usefulPat.pitch = pitchOffset;
	usefulPat.roll = rollOffset;
	usefulPat.heading+=headingOffset;
	if(usefulPat.heading >= 360){
	  usefulPat.heading-=360;
	}
	else if(usefulPat.heading < 0){
	  usefulPat.heading+=360;
	}

	Double_t thetaExpected, phiExpected;
	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected *= TMath::RadToDeg();

	Double_t phi1 = geom->getAntPhiPositionRelToAftFore(ant1)*TMath::RadToDeg();
	Double_t phi2 = geom->getAntPhiPositionRelToAftFore(ant2)*TMath::RadToDeg();
	Double_t deltaAngleDeg1 = RootTools::getDeltaAngleDeg(phi1, phiExpected);
	Double_t deltaAngleDeg2 = RootTools::getDeltaAngleDeg(phi2, phiExpected);

	if(TMath::Abs(deltaAngleDeg1) < maxDeltaAngleDeg && TMath::Abs(deltaAngleDeg2) < maxDeltaAngleDeg){
	
	  Double_t dt_e = usefulPat.getDeltaTExpected(ant2, ant1,
						      AnitaLocations::LONGITUDE_WAIS,
						      AnitaLocations::LATITUDE_WAIS,
						      AnitaLocations::ALTITUDE_WAIS);

	  Double_t diff = dt_e-correlationDeltaTs[combo];
	  sumOfSquaredDiffs += diff*diff;
	  count++;
	}
      }
    }
  }
  sumOfSquaredDiffs/=count;
  sumOfSquaredDiffs = TMath::Sqrt(sumOfSquaredDiffs);
  std::cout << pitchOffset << "\t" << rollOffset << "\t" << headingOffset << "\t" << sumOfSquaredDiffs << std::endl;
  return sumOfSquaredDiffs;

}
