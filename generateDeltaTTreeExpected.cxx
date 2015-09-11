// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Want this program to do three things.
             1.) make deltaTs expected for Prioritizerd antenna positions (i.e. feed positions)
             2.) make deltaTs expected for photogrammetry antenna positions from Linda et. al.
             3.) Somehow find the "best fit" antenna positions given the set of deltaTs for the event.
                 How to do 3?
*************************************************************************************************************** */


////////////////
// Have many, many deltaTs and only a few antenna positions. (But highly correlated information)
// Really only have signal arrival time of each channel... i.e. 48 ts.
// All deltaTs are differences between these numbers...
// 
/////////////


#include <TFile.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>


#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

int main(int argc, char *argv[])
{

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;    
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1] << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;
  const Double_t maxDeltaTriggerTimeNs = 1200;
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2000;
  TString outFileName = TString::Format("%s_run%d-%dPlots.root", argv[0], firstRun, lastRun);
  TFile* outFile = new TFile(outFileName, "recreate");
  TString title1 = TString::Format("Distance from WAIS divide runs %d - %d; run; Distance (km)",
				   firstRun, lastRun);
  TH2D* hDistanceFromWais = new TH2D("hDistanceFromWais", title1,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDistBins, 0, maxDistKm);

  const Double_t minDiff = -2e3;
  const Double_t maxDiff = 2e3;
  const Int_t numDiffBins = 256;
  TString title2 = TString::Format("Event TriggerTimeNs difference from expected WAIS arrival time runs %d - %d; run; #Delta t_{received - expected} (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromWais = new TH2D("hTriggerTimeNs", title2,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDiffBins, minDiff, maxDiff);

  TTree* deltaTTree = new TTree("deltaTTree", "deltaTTree");
  UInt_t eventNumber;
  Double_t heading;
  Double_t thetaExpected;
  Double_t phiExpected;
  UInt_t triggerTimeNs;
  UInt_t triggerTimeNsExpected;    
  std::vector<std::vector<Double_t> >* correlationDeltaTs = NULL;
  std::vector<std::vector<Double_t> >* correlationValues = NULL;  
  std::vector<Double_t> * deltaPhiDeg = NULL;
  deltaTTree->Branch("eventNumber", &eventNumber);
  deltaTTree->Branch("correlationDeltaTs", &correlationDeltaTs);
  deltaTTree->Branch("correlationValues", &correlationValues);  
  deltaTTree->Branch("thetaExpected", &thetaExpected);
  deltaTTree->Branch("phiExpected", &phiExpected);
  deltaTTree->Branch("deltaPhiDeg", &deltaPhiDeg);  
  deltaTTree->Branch("heading", &heading);
  deltaTTree->Branch("triggerTimeNs", &triggerTimeNs);
  deltaTTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);  
  
  Int_t upsampleFactor = 32;
  CrossCorrelator* cc = new CrossCorrelator(upsampleFactor);
  
  // (*thetaExpected) = std::vector<Double_t>(NUM_COMBOS, 0);
  // (*phiExpected) = std::vector<Double_t>(NUM_COMBOS, 0);
  (*deltaPhiDeg) = std::vector<Double_t>(NUM_SEAVEYS, 0);  
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){

	eventNumber = header->eventNumber;
	Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
	hDistanceFromWais->Fill(header->run, distKm);
	hDeltaTFromWais->Fill(header->run, deltaTriggerTimeNs);

	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	usefulPat.getThetaAndPhiWaveAnita3WaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=TMath::RadToDeg();

	for(int ant=0; ant<NUM_SEAVEYS; ant++){
	  deltaPhiDeg->at(ant) = RootTools::getDeltaAngleDeg(phiExpected, cc->phiArrayDeg[ant]);
	}

	cc->correlateEvent(usefulEvent);

	(*correlationDeltaTs) = cc->getMaxCorrelationTimes();
	(*correlationValues) = cc->getMaxCorrelationValues();	

	if(eventNumber==60832576){
	  TGraph* gr = cc->getCrossCorrelationGraph(AnitaPol::kHorizontal, 0, 16);
	  gr->SetName("gr60832576_0_16");
	  gr->Write();
	  delete gr;
	}
	// std::cout << correlationDeltaTs->at(0).at(0) << "\t" << correlationDeltaTs->at(0).size() << "\t" << correlationDeltaTs->size() << std::endl;
	
	delete usefulEvent;

	deltaTTree->Fill();
	
      }
    }
    p++;
  }
  delete cc;
  

  outFile->Write();
  outFile->Close();

  return 0;
}
