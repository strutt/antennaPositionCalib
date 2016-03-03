// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to debug deltaT_expected in CrossCorrelator
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"
#include "FFTtools.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"

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
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;

  const Double_t maxDeltaTriggerTimeNs = 1200;
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  // Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, AnitaPol::kVertical);    
  // if(insertion > 0){
  //   std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
  //   return 1;
  // }

  // insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);
  // if(insertion > 0){
  //   std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
  //   return 1;
  // }

  
  CrossCorrelator* cc = new CrossCorrelator();

  cc->insertPhotogrammetryGeometry();

  const Double_t corDt = cc->correlationDeltaT;
  // cc->kZeroChannel16BH = true;
  // cc->kUseAbbyCombinatorics = 1;
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);  
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  
  lindaFileNameReference->Write();
  
  TTree* corTree = new TTree("corTree", "corTree");
  Double_t phiExpected;
  Double_t thetaExpected;
  corTree->Branch("phiExpected", &phiExpected);
  corTree->Branch("thetaExpected", &thetaExpected);

  Double_t corVals[NUM_COMBOS];
  Double_t corValsClose[NUM_COMBOS];
  Double_t deltaTExpected[NUM_COMBOS];
  Double_t deltaTMeasured[NUM_COMBOS];
  Double_t deltaTMeasured2[NUM_COMBOS];  
  Double_t deltaTMeasuredClose[NUM_COMBOS];
  corTree->Branch(TString::Format("corVals[%d]", NUM_COMBOS), corVals);
  corTree->Branch(TString::Format("corValsClose[%d]", NUM_COMBOS), corValsClose);
  corTree->Branch(TString::Format("deltaTExpected[%d]", NUM_COMBOS), deltaTExpected);
  corTree->Branch(TString::Format("deltaTMeasured[%d]", NUM_COMBOS), deltaTMeasured);
  corTree->Branch(TString::Format("deltaTMeasured2[%d]", NUM_COMBOS), deltaTMeasured2);  
  corTree->Branch(TString::Format("deltaTMeasuredClose[%d]", NUM_COMBOS), deltaTMeasuredClose);

  Double_t heading = 0;
  UInt_t l3TrigPattern = 0;
  UInt_t l3TrigPatternH = 0;
  UInt_t eventNumber = 0;
  UInt_t realTime = 0;
  Int_t run = 0;
  corTree->Branch("heading", &heading);
  corTree->Branch("l3TrigPattern", &l3TrigPattern);
  corTree->Branch("l3TrigPatternH", &l3TrigPatternH);
  corTree->Branch("eventNumber", &eventNumber);
  corTree->Branch("run", &run);
  corTree->Branch("realTime", &realTime);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 1000; //1000; //1000; //10000;
  Long64_t startEntry = 0; //
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    // if(header->eventNumber != 60832108) continue;


    // Int_t numPhiTrigs = 0;
    // for(int phi=0; phi<NUM_PHI; phi++){
    //   Int_t state = RootTools::getBit(phi, l3TrigPatternH);
    //   if(state > 0){
    // 	numPhiTrigs++;
    //   }
    // }

    if((header->trigType & 1)==1){// && numPhiTrigs==2){// && trigIWant > 0){
      gpsChain->GetEntry(entry);
      
      UsefulAdu5Pat usefulPat(pat);
      Double_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      Double_t triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){

	calEventChain->GetEntry(entry);

	l3TrigPattern = header->l3TrigPattern;
	l3TrigPatternH = header->l3TrigPatternH;
	eventNumber = header->eventNumber;
	realTime = header->realTime;
	run = firstRun;
	heading = pat->heading;

	if(pat->heading != usefulPat.heading){
	  std::cerr << "uh oh! " << pat->heading << "\t" << usefulPat.heading << std::endl;
	}
	
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();

	cc->correlateEvent(usefulEvent, pol);
	cc->fillCombosToUseIfNeeded(CrossCorrelator::kTriggered, l3TrigPatternH);
	cc->doUpsampledCrossCorrelationsThreaded(pol, l3TrigPatternH);

	for(int combo=0; combo<NUM_COMBOS; combo++){
	  corVals[combo] = -999;
	  corValsClose[combo] = -999;
	  deltaTExpected[combo] = -999;
	  deltaTMeasured[combo] = -999;
	  deltaTMeasured2[combo] = -999;	  
	  deltaTMeasuredClose[combo] = -999;	  
	}
	
	const std::vector<Int_t>& combos = cc->combosToUseTriggered[l3TrigPatternH];  

	// std::cout << combos.size() << std::endl;
	
	Double_t phiWave = phiExpected*TMath::DegToRad();
	Double_t thetaWave = thetaExpected*TMath::DegToRad();

	for(UInt_t comboInd=0; comboInd < combos.size(); comboInd++){
	  const Int_t combo = combos.at(comboInd);
	  
	  Int_t ant1 = cc->comboToAnt1s.at(combo);
	  Int_t ant2 = cc->comboToAnt2s.at(combo);

	  TGraph* gr1 = cc->grs[pol][ant1];
	  TGraph* gr2 = cc->grs[pol][ant2];

	  TGraph* grCor = FFTtools::getInterpolatedCorrelationGraph(gr2, gr1, corDt);
	  Int_t peakBin = FFTtools::getPeakBin(grCor);
	  deltaTMeasured2[combo] = grCor->GetX()[peakBin];

	  delete grCor;
	  
	  deltaTExpected[combo] = cc->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	  // deltaTExpected[combo] = cc->getDeltaTExpected(pol, ant1, ant2, phiWave, -1*thetaWave);	  
	  Double_t dtExp = deltaTExpected[combo];
	  
	  Double_t tToExp = 1e9;
	  for(int samp=0; samp < cc->numSamplesUpsampled; samp++){
	    Int_t samp0 = samp - 1;
	    Int_t samp1 = samp;
	    Int_t samp2 = samp + 1;

	    samp0 = samp0 < 0 ? samp0 + cc->numSamplesUpsampled : samp0;
	    samp2 = samp2 >= cc->numSamplesUpsampled ? samp2 - cc->numSamplesUpsampled : samp2;
      
	    if(cc->crossCorrelationsUpsampled[pol][combo][samp0] < cc->crossCorrelationsUpsampled[pol][combo][samp1]
	       && cc->crossCorrelationsUpsampled[pol][combo][samp2] < cc->crossCorrelationsUpsampled[pol][combo][samp1]){
	      // Then it's a local maximum
	      // Get time
	      Double_t dtMaxima = samp1 < cc->numSamplesUpsampled/2 ? samp1*corDt : (samp1-cc->numSamplesUpsampled)*corDt;

	      // Is local maximum closest to expected time?
	      if(TMath::Abs(dtMaxima - dtExp) < tToExp){
		tToExp = dtMaxima - dtExp;
		corValsClose[combo] = cc->crossCorrelationsUpsampled[pol][combo][samp1];
		// corInd[combo] = samp1;
		deltaTMeasuredClose[combo] = dtMaxima;
	      }
	    
	      // Is maxima global maximum?
	      if(cc->crossCorrelationsUpsampled[pol][combo][samp1] > corVals[combo]){
		corVals[combo] = cc->crossCorrelationsUpsampled[pol][combo][samp1];
		deltaTMeasured[combo] = dtMaxima;
	      }
	    }
	  }

	}

	delete usefulEvent;
	corTree->Fill();
      }	
    }    
    p++;
  }
  
  outFile->Write();
  outFile->Close();
  return 0;
}


