// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to reconstruct HPol pulses from Wais Divide.
********************************************************************************************************* */

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
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;
  const Double_t maxDeltaTriggerTimeNs = 1200;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  // Photogrammetry positions
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(Int_t surf=0; surf<NUM_SURF; surf++){
    for(Int_t chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
    }
  }
  TString lindaFileName = "photogrammetryNoExtraDelay";
  
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_2015_12_17_time_18_22_28.txt", pol);  
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_LDBHPOL_2015_12_17_time_17_38_21.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_LDBHPOL_2016_01_11_time_11_40_27.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_LDBHPOL_2016_01_11_time_14_42_08.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_2015_10_13_time_14_30_54.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_2016_01_11_time_16_01_54.txt" pol);  
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_NEW3_2016_01_11_time_19_30_31.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_NEW8_2016_01_18_time_17_10_01.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_NEW7_2016_01_18_time_18_18_35.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11.txt", pol);
  // CrossCorrelator::directlyInsertGeometry("newLindaNumbers_LDBHPOL_NEW10_2016_01_19_time_20_38_29.txt", pol);

  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11.txt";
  // TString lindaFileName = "newLindaNumbers_LDBHPOL_NEW10_cosminV2_2016_01_21_time_14_03_25.txt";
  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV2_2016_01_21_time_14_53_41.txt";
  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV3_2016_01_25_time_12_24_25.txt";
  // TString lindaFileName = "newLindaNumbers_LDBHPOL_NEW10_cosminV3_2016_01_25_time_11_33_16.txt";

  // TString lindaFileName = "photogrammetryButWithZeroed16BH";
  // TString lindaFileName = "newBensNumbers_rotatingByHand.txt";
  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  // Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);  
  // if(insertion > 0){
  //   std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
  //   return 1;
  // }

  CrossCorrelator* cc = new CrossCorrelator();  
    
  // cc->kZeroChannel16BH = true;
  
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
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  

  lindaFileNameReference->Write();

  TTree* angResTree = new TTree("angResTree", "angResTree");

  Double_t globalPeak = 0;
  Double_t globalPhiDeg = 0;
  Double_t globalThetaDeg = 0;
  Double_t triggeredPeak = 0;
  Double_t triggeredPhiDeg = 0;
  Double_t triggeredThetaDeg = 0;
  Double_t zoomPeak = 0;
  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t deltaThetaDeg = 0;
  Double_t deltaPhiDeg = 0;
  // Double_t zoomPeak2 = 0;
  // Double_t zoomPhiDeg2 = 0;
  // Double_t zoomThetaDeg2 = 0;
  // Double_t deltaThetaDeg2 = 0;
  // Double_t deltaPhiDeg2 = 0;
  // Double_t zoomPeak3 = 0;
  // Double_t zoomPhiDeg3 = 0;
  // Double_t zoomThetaDeg3 = 0;
  // Double_t deltaThetaDeg3 = 0;
  // Double_t deltaPhiDeg3 = 0;
  // Double_t zoomPeak4 = 0;
  // Double_t zoomPhiDeg4 = 0;
  // Double_t zoomThetaDeg4 = 0;
  // Double_t deltaThetaDeg4 = 0;
  // Double_t deltaPhiDeg4 = 0;
  // Double_t zoomPeak5 = 0;
  // Double_t zoomPhiDeg5 = 0;
  // Double_t zoomThetaDeg5 = 0;
  // Double_t deltaThetaDeg5 = 0;
  // Double_t deltaPhiDeg5 = 0;
  // Double_t zoomPeak6 = 0;
  // Double_t zoomPhiDeg6 = 0;
  // Double_t zoomThetaDeg6 = 0;
  // Double_t deltaThetaDeg6 = 0;
  // Double_t deltaPhiDeg6 = 0;
  // Double_t zoomPeak7 = 0;
  // Double_t zoomPhiDeg7 = 0;
  // Double_t zoomThetaDeg7 = 0;
  // Double_t deltaThetaDeg7 = 0;
  // Double_t deltaPhiDeg7 = 0;
  // Double_t zoomPeak8 = 0;
  // Double_t zoomPhiDeg8 = 0;
  // Double_t zoomThetaDeg8 = 0;
  // Double_t deltaThetaDeg8 = 0;
  // Double_t deltaPhiDeg8 = 0;

  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
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
  
  angResTree->Branch("globalPeak", &globalPeak);
  angResTree->Branch("globalPhiDeg", &globalPhiDeg);
  angResTree->Branch("globalThetaDeg", &globalThetaDeg);

  angResTree->Branch("triggeredPeak", &triggeredPeak);
  angResTree->Branch("triggeredPhiDeg", &triggeredPhiDeg);
  angResTree->Branch("triggeredThetaDeg", &triggeredThetaDeg);

  angResTree->Branch("zoomPeak", &zoomPeak);
  angResTree->Branch("zoomPhiDeg", &zoomPhiDeg);
  angResTree->Branch("zoomThetaDeg", &zoomThetaDeg);
  angResTree->Branch("deltaPhiDeg", &deltaPhiDeg);
  angResTree->Branch("deltaThetaDeg", &deltaThetaDeg);

  // angResTree->Branch("zoomPeak2", &zoomPeak2);
  // angResTree->Branch("zoomPhiDeg2", &zoomPhiDeg2);
  // angResTree->Branch("zoomThetaDeg2", &zoomThetaDeg2);
  // angResTree->Branch("deltaPhiDeg2", &deltaPhiDeg2);
  // angResTree->Branch("deltaThetaDeg2", &deltaThetaDeg2);

  // angResTree->Branch("zoomPeak3", &zoomPeak3);
  // angResTree->Branch("zoomPhiDeg3", &zoomPhiDeg3);
  // angResTree->Branch("zoomThetaDeg3", &zoomThetaDeg3);
  // angResTree->Branch("deltaPhiDeg3", &deltaPhiDeg3);
  // angResTree->Branch("deltaThetaDeg3", &deltaThetaDeg3);

  // angResTree->Branch("zoomPeak4", &zoomPeak4);
  // angResTree->Branch("zoomPhiDeg4", &zoomPhiDeg4);
  // angResTree->Branch("zoomThetaDeg4", &zoomThetaDeg4);
  // angResTree->Branch("deltaPhiDeg4", &deltaPhiDeg4);
  // angResTree->Branch("deltaThetaDeg4", &deltaThetaDeg4);

  // angResTree->Branch("zoomPeak5", &zoomPeak5);
  // angResTree->Branch("zoomPhiDeg5", &zoomPhiDeg5);
  // angResTree->Branch("zoomThetaDeg5", &zoomThetaDeg5);
  // angResTree->Branch("deltaPhiDeg5", &deltaPhiDeg5);
  // angResTree->Branch("deltaThetaDeg5", &deltaThetaDeg5);

  // angResTree->Branch("zoomPeak6", &zoomPeak6);
  // angResTree->Branch("zoomPhiDeg6", &zoomPhiDeg6);
  // angResTree->Branch("zoomThetaDeg6", &zoomThetaDeg6);
  // angResTree->Branch("deltaPhiDeg6", &deltaPhiDeg6);
  // angResTree->Branch("deltaThetaDeg6", &deltaThetaDeg6);

  // angResTree->Branch("zoomPeak7", &zoomPeak7);
  // angResTree->Branch("zoomPhiDeg7", &zoomPhiDeg7);
  // angResTree->Branch("zoomThetaDeg7", &zoomThetaDeg7);
  // angResTree->Branch("deltaPhiDeg7", &deltaPhiDeg7);
  // angResTree->Branch("deltaThetaDeg7", &deltaThetaDeg7);

  // angResTree->Branch("zoomPeak8", &zoomPeak8);
  // angResTree->Branch("zoomPhiDeg8", &zoomPhiDeg8);
  // angResTree->Branch("zoomThetaDeg8", &zoomThetaDeg8);
  // angResTree->Branch("deltaPhiDeg8", &deltaPhiDeg8);
  // angResTree->Branch("deltaThetaDeg8", &deltaThetaDeg8);
  
  angResTree->Branch("phiSectorOfPeak", &phiSectorOfPeak);
  angResTree->Branch("hackyL3Trig", &hackyL3Trig);
  
  angResTree->Branch("thetaExpected", &thetaExpected);
  angResTree->Branch("phiExpected", &phiExpected);
  angResTree->Branch("triggerTimeNs", &triggerTimeNs);
  angResTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  angResTree->Branch("heading", &heading);
  angResTree->Branch("l3TrigPattern", &l3TrigPattern);
  angResTree->Branch("l3TrigPatternH", &l3TrigPatternH);
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);

  angResTree->Branch("mrms", &mrms);
  angResTree->Branch("brms", &brms);
  angResTree->Branch("realTime", &realTime);

  Double_t maxDPhiDeg;
  angResTree->Branch("maxDPhiDeg", &maxDPhiDeg);
  // Double_t maxDPhiDeg2;
  // angResTree->Branch("maxDPhiDeg2", &maxDPhiDeg2);
  // Double_t maxDPhiDeg3;
  // angResTree->Branch("maxDPhiDeg3", &maxDPhiDeg3);
  // Double_t maxDPhiDeg4;
  // angResTree->Branch("maxDPhiDeg4", &maxDPhiDeg4);
  // Double_t maxDPhiDeg5;
  // angResTree->Branch("maxDPhiDeg5", &maxDPhiDeg5);
  // Double_t maxDPhiDeg6;
  // angResTree->Branch("maxDPhiDeg6", &maxDPhiDeg6);
  // Double_t maxDPhiDeg7;
  // angResTree->Branch("maxDPhiDeg7", &maxDPhiDeg7);
  // Double_t maxDPhiDeg8;
  // angResTree->Branch("maxDPhiDeg8", &maxDPhiDeg8);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //10000;
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
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
      // if(true){
	eventNumber = header->eventNumber;
	run = header->run;
	
	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	mrms = usefulPat.mrms;
	brms = usefulPat.brms;
	realTime = header->realTime;
	l3TrigPattern = header->l3TrigPattern;
	l3TrigPatternH = header->l3TrigPatternH;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();
	
	phiExpected = phiExpected < 0 ? phiExpected + 360 : phiExpected;
	phiExpected = phiExpected >= 360 ? phiExpected - 360 : phiExpected;

	cc->correlateEvent(usefulEvent, pol);
	
	TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);
	hGlobalImageH->SetName(TString::Format("hGlobalImageH_%u", eventNumber));
	// globalPhiDeg = -999;
	// globalThetaDeg = -999;
	
	// generic triggered, for getting close to target!

	// cc->maxDPhiDeg = 0;
	
	// cc->kUseOffAxisDelay = 0;
	// cc->kDeltaPhiSect = 2;
	// TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
	// 					       triggeredThetaDeg, l3TrigPatternH);

	// triggeredPhiDeg = triggeredPhiDeg < 0 ? triggeredPhiDeg + 360 : triggeredPhiDeg;
	// triggeredPhiDeg = triggeredPhiDeg >= 360 ? triggeredPhiDeg - 360 : triggeredPhiDeg;

	phiSectorOfPeak = -1;
	Double_t bestDeltaPhiOfPeakToAnt = 360;
	for(int ant=0; ant < NUM_SEAVEYS; ant++){
	  Double_t phiOfAnt = cc->phiArrayDeg[pol].at(ant);
	  Double_t deltaPhiOfPeakToAnt = TMath::Abs(RootTools::getDeltaAngleDeg(phiOfAnt, triggeredPhiDeg));
	  if(deltaPhiOfPeakToAnt < bestDeltaPhiOfPeakToAnt){
	    bestDeltaPhiOfPeakToAnt = deltaPhiOfPeakToAnt;
	    phiSectorOfPeak = (ant % NUM_PHI);
	  }
	}
	hackyL3Trig = (1 << phiSectorOfPeak);

	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						  zoomThetaDeg,
						  globalPhiDeg, globalThetaDeg);
						  // triggeredPhiDeg, triggeredThetaDeg);

	zoomPhiDeg = zoomPhiDeg < 0 ? zoomPhiDeg + 360 : zoomPhiDeg;
	zoomPhiDeg = zoomPhiDeg >= 360 ? zoomPhiDeg - 360 : zoomPhiDeg;

	deltaPhiDeg = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg);
	deltaThetaDeg = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg);

	TString name1 = TString::Format("hTriggeredZoomImageH%u_1", eventNumber);
	hZoomedImageH->SetName(name1);

	maxDPhiDeg = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kDeltaPhiSect = 1;

	// TH2D* hZoomedImageH2 = cc->makeZoomedImage(pol, zoomPeak2, zoomPhiDeg2,
	// 					 zoomThetaDeg2, triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg2 = zoomPhiDeg2 < 0 ? zoomPhiDeg2 + 360 : zoomPhiDeg2;
	// zoomPhiDeg2 = zoomPhiDeg2 >= 360 ? zoomPhiDeg2 - 360 : zoomPhiDeg2;
	
	// deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg2);
	// deltaThetaDeg2 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg2);

	// TString name2 = TString::Format("hTriggeredZoomImageH%u_2", eventNumber);
	// hZoomedImageH2->SetName(name2);


	// maxDPhiDeg2 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kDeltaPhiSect = 2;
	// cc->kUseOffAxisDelay = 1;	

	// TH2D* hZoomedImageH3 = cc->makeZoomedImage(pol, zoomPeak3, zoomPhiDeg3,
	// 					 zoomThetaDeg3, hackyL3Trig,
	// 					 triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg3 = zoomPhiDeg3 < 0 ? zoomPhiDeg3 + 360 : zoomPhiDeg3;
	// zoomPhiDeg3 = zoomPhiDeg3 >= 360 ? zoomPhiDeg3 - 360 : zoomPhiDeg3;
	
	// deltaPhiDeg3 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg3);
	// deltaThetaDeg3 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg3);

	// TString name3 = TString::Format("hTriggeredZoomImageH%u_3", eventNumber);
	// hZoomedImageH3->SetName(name3);

	

	// maxDPhiDeg3 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kDeltaPhiSect = 1;

	// TH2D* hZoomedImageH4 = cc->makeZoomedImage(pol, zoomPeak4, zoomPhiDeg4,
	// 					 zoomThetaDeg4, hackyL3Trig,
	// 					 triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg4 = zoomPhiDeg4 < 0 ? zoomPhiDeg4 + 360 : zoomPhiDeg4;
	// zoomPhiDeg4 = zoomPhiDeg4 >= 360 ? zoomPhiDeg4 - 360 : zoomPhiDeg4;
	
	// deltaPhiDeg4 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg4);
	// deltaThetaDeg4 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg4);


	// TString name4 = TString::Format("hTriggeredZoomImageH%u_4", eventNumber);
	// hZoomedImageH4->SetName(name4);


	// maxDPhiDeg4 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;	
	// cc->kUseOffAxisDelay = 0;
	// cc->kDeltaPhiSect = -2;

	// TH2D* hZoomedImageH5 = cc->makeZoomedImage(pol, zoomPeak5, zoomPhiDeg5,
	// 					 zoomThetaDeg5, triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg5 = zoomPhiDeg5 < 0 ? zoomPhiDeg5 + 360 : zoomPhiDeg5;
	// zoomPhiDeg5 = zoomPhiDeg5 >= 360 ? zoomPhiDeg5 - 360 : zoomPhiDeg5;
	
	// deltaPhiDeg5 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg5);
	// deltaThetaDeg5 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg5);	


	// TString name5 = TString::Format("hTriggeredZoomImageH%u_5", eventNumber);
	// hZoomedImageH5->SetName(name5);


	
	
	// maxDPhiDeg5 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kDeltaPhiSect = -1;

	// TH2D* hZoomedImageH6 = cc->makeZoomedImage(pol, zoomPeak6, zoomPhiDeg6,
	// 					 zoomThetaDeg6, triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg6 = zoomPhiDeg6 < 0 ? zoomPhiDeg6 + 360 : zoomPhiDeg6;
	// zoomPhiDeg6 = zoomPhiDeg6 >= 360 ? zoomPhiDeg6 - 360 : zoomPhiDeg6;
	
	// deltaPhiDeg6 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg6);
	// deltaThetaDeg6 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg6);	


	// TString name6 = TString::Format("hTriggeredZoomImageH%u_6", eventNumber);
	// hZoomedImageH6->SetName(name6);



	// maxDPhiDeg6 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kUseOffAxisDelay = 1;
	// cc->kDeltaPhiSect = -2;

	// TH2D* hZoomedImageH7 = cc->makeZoomedImage(pol, zoomPeak7, zoomPhiDeg7,
	// 					 zoomThetaDeg7, hackyL3Trig,
	// 					 triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg7 = zoomPhiDeg7 < 0 ? zoomPhiDeg7 + 360 : zoomPhiDeg7;
	// zoomPhiDeg7 = zoomPhiDeg7 >= 360 ? zoomPhiDeg7 - 360 : zoomPhiDeg7;
	
	// deltaPhiDeg7 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg7);
	// deltaThetaDeg7 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg7);	


	// TString name7 = TString::Format("hTriggeredZoomImageH%u_7", eventNumber);
	// hZoomedImageH7->SetName(name7);


	
	
	// maxDPhiDeg7 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	// cc->kDeltaPhiSect = -1;

	// TH2D* hZoomedImageH8 = cc->makeZoomedImage(pol, zoomPeak8, zoomPhiDeg8,
	// 					 zoomThetaDeg8, hackyL3Trig,
	// 					 triggeredPhiDeg, triggeredThetaDeg);

	// zoomPhiDeg8 = zoomPhiDeg8 < 0 ? zoomPhiDeg8 + 360 : zoomPhiDeg8;
	// zoomPhiDeg8 = zoomPhiDeg8 >= 360 ? zoomPhiDeg8 - 360 : zoomPhiDeg8;
	
	// deltaPhiDeg8 = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg8);
	// deltaThetaDeg8 = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg8);	


	// TString name8 = TString::Format("hTriggeredZoomImageH%u_8", eventNumber);
	// hZoomedImageH8->SetName(name8);



	// maxDPhiDeg8 = cc->maxDPhiDeg;
	// cc->maxDPhiDeg = 0;
	
	if(numSaved < maxToSave){
	  numSaved++;
	}
	else{
	  // delete hTriggeredImageH;
	  delete hGlobalImageH;
	  delete hZoomedImageH;
	  // delete hZoomedImageH2;
	  // delete hZoomedImageH3;
	  // delete hZoomedImageH4;
	  // delete hZoomedImageH5;
	  // delete hZoomedImageH6;	  
	  // delete hZoomedImageH7;
	  // delete hZoomedImageH8;	  
	}

	angResTree->Fill();

	delete usefulEvent;
      }
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
