// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to reconstruct VPol pulses from LDB.
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
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;
  // const Int_t cutTimeNs = Int_t(60e3);
  const Int_t cutTimeNs = Int_t(2e3);  
  
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

  //  Photogrammetry positions
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(Int_t surf=0; surf<NUM_SURF; surf++){
    for(Int_t chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
    }
  }
  TString lindaFileName = "photogrammetryNoExtraDelay";

  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_LDBHPOL_2015_12_17_time_17_38_21.txt",
  // 					  pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_VPOL_10kVSeavey_2015_12_17_time_17_41_28.txt",  pol);  
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_VPOL_10kVSeavey_TEEEEEEEST_2016_01_11_time_14_56_23.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_2015_10_13_time_14_30_54.txt", pol); 
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_2016_01_19_time_20_39_33", pol);

  // TString lindaFileName = "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_2016_01_19_time_20_39_33.txt";
  // TString lindaFileName = "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_cosminV2_2016_01_21_time_14_30_18.txt";
  // TString lindaFileName = "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_cosminV3_2016_01_25_time_11_44_07.txt";

  // TString lindaFileName = "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW12_bugFixed_reFit_2016_02_18_time_18_38_23.txt";

  // CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);
  
  CrossCorrelator* cc = new CrossCorrelator();
  
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%d.root", run, run);    
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);    
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);  
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  // gpsChain->BuildIndex("realTime");  
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);
  
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  lindaFileNameReference->Write();
  
  const Double_t sourceLat = - (77 + (51.23017/60)); // OLD NUMBERS
  // const Double_t sourceLat = - (77 + 51./60 + 44.16/3600); // From Google Earth 77°51'44.16"S
  const Double_t sourceLon = +(167 + (12.16908/60));// OLD NUMBERS
  // const Double_t sourceLon = +(167 + 2./60 + 28.83/3600); // From Google Earth 167° 2'28.83"E
  const Double_t sourceAlt = 0;
  
  TTree* angResTree = new TTree("angResTree", "angResTree");

  Double_t globalPeak = 0;
  Double_t globalPhiDeg = 0;
  Double_t globalThetaDeg = 0;
  // Double_t triggeredPeak = 0;
  // Double_t triggeredPhiDeg = 0;
  // Double_t triggeredThetaDeg = 0;
  Double_t zoomPeak = 0;
  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  Double_t deltaThetaDeg = 0;
  Double_t deltaPhiDeg = 0;
  UInt_t eventNumber = 0;
  UInt_t triggerTimeNs = 0;
  UInt_t triggerTimeNsExpected = 0;
  Double_t heading = 0;
  UInt_t l3TrigPattern = 0;
  UInt_t l3TrigPatternH = 0;  
  Int_t run = 0;  
  Double_t mrms = 0;
  Double_t brms = 0;
  UInt_t realTime = 0;    
  Int_t phiSectorOfPeak = 0;
  UShort_t hackyL3Trig = 0;

  // std::vector<Double_t>* deltaPhiDeg = NULL;

  angResTree->Branch("globalPeak", &globalPeak);
  angResTree->Branch("globalPhiDeg", &globalPhiDeg);
  angResTree->Branch("globalThetaDeg", &globalThetaDeg);

  // angResTree->Branch("triggeredPeak", &triggeredPeak);
  // angResTree->Branch("triggeredPhiDeg", &triggeredPhiDeg);
  // angResTree->Branch("triggeredThetaDeg", &triggeredThetaDeg);

  angResTree->Branch("zoomPeak", &zoomPeak);
  angResTree->Branch("zoomPhiDeg", &zoomPhiDeg);
  angResTree->Branch("zoomThetaDeg", &zoomThetaDeg);

  angResTree->Branch("deltaPhiDeg", &deltaPhiDeg);
  angResTree->Branch("deltaThetaDeg", &deltaThetaDeg);

  angResTree->Branch("thetaExpected", &thetaExpected);
  angResTree->Branch("phiExpected", &phiExpected);
  angResTree->Branch("triggerTimeNs", &triggerTimeNs);
  angResTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  angResTree->Branch("heading", &heading);
  angResTree->Branch("l3TrigPattern", &l3TrigPattern);
  angResTree->Branch("l3TrigPatternH", &l3TrigPatternH);  
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);

  angResTree->Branch("phiSectorOfPeak", &phiSectorOfPeak);
  angResTree->Branch("hackyL3Trig", &hackyL3Trig);
  
  angResTree->Branch("mrms", &mrms);
  angResTree->Branch("brms", &brms);
  angResTree->Branch("realTime", &realTime);  

  UInt_t timedelay[10];
  Int_t shortdelay = 3150;
  int ndelays = -1;
  Int_t deltaTns_corr[5];
  if (pol==AnitaPol::kVertical){
    timedelay[0] = Int_t(25e6);
    timedelay[1] = Int_t(225e6);
    timedelay[2] = Int_t(425e6);
    timedelay[3] = Int_t(625e6);
    timedelay[4] = Int_t(825e6); // in ms
    if (firstRun>153){
      timedelay[0] = Int_t(50e6);
      timedelay[1] = Int_t(250e6);
      timedelay[2] = Int_t(450e6);
      timedelay[3] = Int_t(650e6);
      timedelay[4] = Int_t(850e6);
    }
    ndelays = 5;
  } 

  Long64_t timedelay_corr[10];
  for (int i=0;i<ndelays;++i){
    timedelay_corr[i]= timedelay[i] - i*shortdelay;
  }


  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0;//5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 5;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    // gpsChain->GetEntryWithIndex(header->realTime);
    gpsChain->GetEntry(entry);
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      triggerTimeNs = header->triggerTimeNs;

      Long64_t minDeltaT = Long64_t(1e11);
      for (int i=0;i<ndelays;++i){
	deltaTns_corr[i]= triggerTimeNsExpected+timedelay_corr[i];
	if (deltaTns_corr[i]>Int_t(1e9)) deltaTns_corr[i]-=Int_t(1e9);
	deltaTns_corr[i] = deltaTns_corr[i] - triggerTimeNs;
	if ( TMath::Abs(deltaTns_corr[i]) < TMath::Abs(minDeltaT) ){
	  minDeltaT = deltaTns_corr[i];
	}
      }
      
     if(TMath::Abs(minDeltaT) < cutTimeNs){
	eventNumber = header->eventNumber;
	run = header->run;

	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	mrms = usefulPat.mrms;
	brms = usefulPat.brms;
	realTime = header->realTime;
	l3TrigPatternH = header->l3TrigPatternH;
	l3TrigPattern = header->l3TrigPattern;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
	globalPhiDeg = usefulEvent->eventNumber;
	
	usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();

	// cc->correlateEvent(usefulEvent, pol);

	cc->getNormalizedInterpolatedTGraphs(usefulEvent, pol);

	// cc->doFFTs(pol);
	for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	  FancyFFTs::doFFT(cc->numSamples, cc->grsResampled[pol][ant]->GetY(), cc->ffts[pol][ant]);

	  // if(cc->kDoSimpleSatelliteFiltering > 0){
	  // cc->simple260MHzSatelliteNotch(pol, ant);
	  // cc->simple370MHzSatelliteNotch(pol, ant);
	  const int numFreqs = FancyFFTs::getNumFreqs(cc->numSamples);
	  const Double_t deltaF_MHz = 1e3/(cc->numSamples*cc->nominalSamplingDeltaT);

	  const Double_t notchLowEdgeMHz = 0;
	  const Double_t notchHighEdgeMHz = 500;

	  const Double_t notchLowEdgeMHz2 = 850;
	  const Double_t notchHighEdgeMHz2 = 1400;

	  for(int freqInd=0; freqInd < numFreqs; freqInd++){
	    // if(deltaF_MHz*freqInd >= notchLowEdgeMHz && deltaF_MHz*freqInd < notchHighEdgeMHz){
	    if((deltaF_MHz*freqInd >= notchLowEdgeMHz && deltaF_MHz*freqInd < notchHighEdgeMHz) || (deltaF_MHz*freqInd >= notchLowEdgeMHz2 && deltaF_MHz*freqInd < notchHighEdgeMHz2)){
	      
 	      // std::cout << pol << "\t" << ant << "\t" << deltaF_MHz << "\t" << deltaF_MHz*freqInd << "\t" << notchLowEdgeMHz << "\t" << notchHighEdgeMHz << "\t" << std::endl;
	      cc->ffts[pol][ant][freqInd].real(0);
	      cc->ffts[pol][ant][freqInd].imag(0);
	    }    
	  }
	    
	  cc->renormalizeFourierDomain(pol, ant);

	  FancyFFTs::zeroPadFFT(cc->ffts[pol][ant],
				cc->fftsPadded[pol][ant],
				cc->numSamples,
				cc->numSamplesUpsampled);
	  
	}	  
	cc->doAllCrossCorrelationsThreaded(pol);
	cc->eventNumber[pol] = usefulEvent->eventNumber;	

	// reconstruct
	cc->reconstruct(pol, cc->coarseMapPeakValues[pol][0],
			cc->coarseMapPeakPhiDegs[pol][0],
			cc->coarseMapPeakThetaDegs[pol][0]);


	TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);

	// TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
	// 					       triggeredThetaDeg, l3TrigPattern);

	phiSectorOfPeak = -1;
	Double_t bestDeltaPhiOfPeakToAnt = 360;
	for(int ant=0; ant < NUM_SEAVEYS; ant++){
	  Double_t phiOfAnt = cc->phiArrayDeg[pol].at(ant);
	  Double_t deltaPhiOfPeakToAnt = TMath::Abs(RootTools::getDeltaAngleDeg(phiOfAnt, globalPhiDeg));
	  // Double_t deltaPhiOfPeakToAnt = TMath::Abs(RootTools::getDeltaAngleDeg(phiOfAnt, triggeredPhiDeg));	  
	  
	  if(deltaPhiOfPeakToAnt < bestDeltaPhiOfPeakToAnt){
	    bestDeltaPhiOfPeakToAnt = deltaPhiOfPeakToAnt;
	    phiSectorOfPeak = (ant % NUM_PHI);
	  }
	}
	hackyL3Trig = (1 << phiSectorOfPeak);
	

	// TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
	// 					 zoomThetaDeg, l3TrigPattern,
	// 					 triggeredPhiDeg, triggeredThetaDeg);
	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						 zoomThetaDeg, hackyL3Trig,
						 // triggeredPhiDeg, triggeredThetaDeg);
						 globalPhiDeg, globalThetaDeg);		
	
	globalPhiDeg = globalPhiDeg < 0 ? globalPhiDeg + 360 : globalPhiDeg;
	globalPhiDeg = globalPhiDeg >= 360 ? globalPhiDeg - 360 : globalPhiDeg;

	// triggeredPhiDeg = triggeredPhiDeg < 0 ? triggeredPhiDeg + 360 : triggeredPhiDeg;
	// triggeredPhiDeg = triggeredPhiDeg >= 360 ? triggeredPhiDeg - 360 : triggeredPhiDeg;

	zoomPhiDeg = zoomPhiDeg < 0 ? zoomPhiDeg + 360 : zoomPhiDeg;
	zoomPhiDeg = zoomPhiDeg >= 360 ? zoomPhiDeg - 360 : zoomPhiDeg;

	phiExpected = phiExpected < 0 ? phiExpected + 360 : phiExpected;
	phiExpected = phiExpected >= 360 ? phiExpected - 360 : phiExpected;

	deltaPhiDeg = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg);
	deltaThetaDeg = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg);
	
	
	if(numSaved < maxToSave){
	  numSaved++;
	}
	else{
	  // delete hTriggeredImageH;
	  delete hGlobalImageH;
	  delete hZoomedImageH;	  
	}

	angResTree->Fill();

	delete usefulEvent;
	
      }
    }
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
