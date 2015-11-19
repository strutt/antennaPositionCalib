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


#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>
#include <AnitaEventCalibrator.h>
#include <AnitaGeomTool.h>

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
  // const Double_t maxDeltaTriggerTimeNs = 1200;

  const Double_t minDeltaTriggerTimeNs = 24.994e6;
  const Double_t maxDeltaTriggerTimeNs = 25.003e6;

  const Double_t sourceLat = - (77 + (51.23017/60));
  const Double_t sourceLon = +(167 + (12.16908/60));
  const Double_t sourceAlt = 0;
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnitaIIINumbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(int surf=0; surf<NUM_SURF; surf++){
    for(int chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
    }
  }

  

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");


  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);  
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");  
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2000;  
  TString outFileName = TString::Format("%s_run%d-%dPlots.root", argv[0], firstRun, lastRun);
  TFile* outFile = new TFile(outFileName, "recreate");    
  TString title1 = TString::Format("Distance from LDB divide runs %d - %d; run; Distance (km)", firstRun, lastRun);
  TH2D* hDistanceFromLdb = new TH2D("hDistanceFromLdb", title1,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDistBins, 0, maxDistKm);


  const Double_t minDiff = -2e3;
  const Double_t maxDiff = 2e3;
  const Int_t numDiffBins = 256;
  TString title2 = TString::Format("Event TriggerTimeNs difference from expected LDB arrival time runs %d - %d; run; #Delta t_{received - expected} (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromLdb = new TH2D("hTriggerTimeNs", title2,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDiffBins, minDiff, maxDiff);

  TTree* deltaTTree = new TTree("deltaTTree", "deltaTTree");
  UInt_t eventNumber;
  Double_t heading;
  Double_t thetaExpected;
  Double_t phiExpected;
  UInt_t triggerTimeNs;
  UInt_t triggerTimeNsExpected;    
  // std::vector<std::vector<Double_t> >* correlationDeltaTs = NULL;
  // std::vector<std::vector<Double_t> >* correlationValues = NULL;  
  // std::vector<Double_t> * deltaPhiDeg = NULL;
  Double_t correlationDeltaTs[NUM_COMBOS] = {0};
  Double_t correlationValues[NUM_COMBOS] = {0};
  Double_t correlationDeltaTsClose[NUM_COMBOS] = {0};
  Double_t correlationValuesClose[NUM_COMBOS] = {0};
  std::vector<Double_t> * deltaPhiDeg = NULL;

  deltaTTree->Branch("eventNumber", &eventNumber);
  deltaTTree->Branch(TString::Format("correlationDeltaTs[%d]", NUM_COMBOS), correlationDeltaTs);
  deltaTTree->Branch(TString::Format("correlationValues[%d]", NUM_COMBOS), correlationValues);
  deltaTTree->Branch(TString::Format("correlationDeltaTsClose[%d]", NUM_COMBOS), correlationDeltaTsClose);
  deltaTTree->Branch(TString::Format("correlationValuesClose[%d]", NUM_COMBOS), correlationValuesClose);
  
  // deltaTTree->Branch("correlationDeltaTs", &correlationDeltaTs,
  // 		     TString::Format("correlationDeltaTs[%d][%d]/D", NUM_POL, NUM_COMBOS));
  // deltaTTree->Branch("correlationValues", &correlationValues,
  // 		     TString::Format("correlationValues[%d][%d]/D", NUM_POL, NUM_COMBOS));
  deltaTTree->Branch("thetaExpected", &thetaExpected);
  deltaTTree->Branch("phiExpected", &phiExpected);
  deltaTTree->Branch("deltaPhiDeg", &deltaPhiDeg);
  deltaTTree->Branch("heading", &heading);
  deltaTTree->Branch("triggerTimeNs", &triggerTimeNs);
  deltaTTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  
  CrossCorrelator* cc = new CrossCorrelator();
  
  // (*thetaExpected) = std::vector<Double_t>(NUM_COMBOS, 0);
  // (*phiExpected) = std::vector<Double_t>(NUM_COMBOS, 0);
  (*deltaPhiDeg) = std::vector<Double_t>(NUM_SEAVEYS, 0);

  const Double_t deltaTSearchLimit = 1; //ns
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      triggerTimeNs = header->triggerTimeNs;

      Int_t a = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected) - 0.0005e6; // + 0.003e8;
      a = a%Int_t(1999969e2);
      
      if(a > 67.5e6){
      	a -= 50.0012e6;
      }
      else if(a > 27.5e6){
      	a -= 25e6;
      }
      if(a > minDeltaTriggerTimeNs && a < maxDeltaTriggerTimeNs){
      
      // UsefulAdu5Pat usefulPat(pat);
      // triggerTimeNsExpected = usefulPat.getLdbDivideTriggerTimeNs();
      // triggerTimeNs = header->triggerTimeNs;
      // Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      // if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
	eventNumber = header->eventNumber;

	// if(eventNumber!=60834309) continue;
	
	Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
	hDistanceFromLdb->Fill(header->run, distKm);
	hDeltaTFromLdb->Fill(header->run, Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected));

	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	// usefulPat.getThetaAndPhiWaveLdbDivide(thetaExpected, phiExpected);
	// usefulPat.getThetaAndPhiWaveLdbDivide(thetaExpected, phiExpected);
	usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaExpected, phiExpected);	

	// cc->correlateEvent(usefulEvent, AnitaPol::kHorizontal);

	// Interface in flux, let's hack our way through this...
	cc->getNormalizedInterpolatedTGraphs(usefulEvent, AnitaPol::kVertical);  
	cc->doFFTs(AnitaPol::kVertical);
	cc->fillCombosToUseIfNeeded(CrossCorrelator::kTriggered, 0xffff);
	cc->doUpsampledCrossCorrelationsThreaded(AnitaPol::kVertical, 0xffff);


	// for(Int_t combo=0; combo<NUM_COMBOS; combo++){	
	//   std::cout << combo << "\t" << cc->crossCorrelationsUpsampled[AnitaPol::kVertical][combo] << std::endl;
	// }
	
	for(Int_t combo=0; combo<NUM_COMBOS; combo++){
	  cc->getMaxUpsampledCorrelationTimeValue(AnitaPol::kVertical, combo,
						  correlationDeltaTs[combo],
						  correlationValues[combo]);
	}
	
	for(Int_t combo=0; combo<NUM_COMBOS; combo++){
	  Int_t ant1 = cc->comboToAnt1s.at(combo);
	  Int_t ant2 = cc->comboToAnt2s.at(combo);
	  // Double_t dtExpected = usefulPat.getDeltaTExpected(ant1, ant2, phiExpected, thetaExpected);
	  Double_t dtExpected = usefulPat.getDeltaTExpected(ant2, ant1, phiExpected, thetaExpected);	  

	  TGraph* gr = cc->getUpsampledCrossCorrelationGraph(AnitaPol::kVertical, ant1, ant2);
	    
	  Double_t minY;
	  Double_t minX;
	  RootTools::getMaxMinWithinLimits(gr, correlationValuesClose[combo],
					   correlationDeltaTsClose[combo],
					   minY, minX,
					   dtExpected-deltaTSearchLimit,
					   dtExpected+deltaTSearchLimit);
	  // std::cout << ant1 << "\t" << ant2 << "\t" << dtExpected << "\t" << correlationDeltaTs[combo] << "\t" << correlationValues[combo] << "\t" << dtExpected-deltaTSearchLimit << "\t" << dtExpected+deltaTSearchLimit << std::endl;

	  
	  if(eventNumber==60834309){
	    gr->SetName(TString::Format("gr%u_%d_%d", eventNumber, ant1, ant2));
	    gr->Write();
	  }
	  
	  delete gr;
	}

	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();
	for(int ant=0; ant<NUM_SEAVEYS; ant++){
	  Double_t antPhiDeg = geom->getAntPhiPositionRelToAftFore(ant)*TMath::RadToDeg();
	  deltaPhiDeg->at(ant) = RootTools::getDeltaAngleDeg(phiExpected, antPhiDeg);
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
