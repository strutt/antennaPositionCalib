// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to find peak cross correlation offsets between antenna pairs in pulses from Wais Divide.
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
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;
  const Double_t maxDeltaTriggerTimeNs = 1200;

  // Double_t phiBefore[NUM_SEAVEYS];
  // Double_t phiAfter[NUM_SEAVEYS];  
  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // for(int ant=0; ant<NUM_SEAVEYS; ant++){
  //   phiBefore[ant] = geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::kVertical);
  //   // phiBefore[ant] += TMath::DegToRad()*45;
  // }

  // Validate Linda's Numbers
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  // TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_2015_12_06_time_18_23_32"; 
  TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_2015_12_17_time_18_22_28";

  TString lindaFileNameTxt = lindaFileName + ".txt";
  Int_t geomVal = CrossCorrelator::directlyInsertGeometry(lindaFileNameTxt, pol);
  
  if(geomVal!=0){
    std::cerr << "Could not find file " << lindaFileNameTxt.Data() << std::endl;
    return 1;
  }
  geomVal = CrossCorrelator::validateGeometry(lindaFileNameTxt, pol);
  if(geomVal!=0){
    std::cerr << "Could not validate Linda's numbers from " << lindaFileNameTxt.Data() << std::endl;
    return 1;
  }

  // for(int ant=0; ant<NUM_SEAVEYS; ant++){
  //   phiAfter[ant] = geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::kVertical);
  //   // phiAfter[ant] += TMath::DegToRad()*45;    

  //   if(ant==47){
  //     std::cout << phiBefore[ant] << "\t" << phiAfter[ant] << "\t" << phiBefore[ant] - phiAfter[ant] << std::endl;
  //     std::cout << "geom = " << geom << std::endl;
  //   }
  // }

  

  CrossCorrelator* cc = new CrossCorrelator();  
  // std::cout << "in cc round 2, " << cc->phiArrayDeg[AnitaPol::kVertical][47]*TMath::DegToRad() << std::endl;
  // return 0;

  
  // Get input

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
  
  OutputConvention oc(argc, argv);
  oc.setSubdirectory(lindaFileName);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

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
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  Double_t deltaThetaDeg = 0;
  Double_t deltaPhiDeg = 0;
  UInt_t eventNumber = 0;
  UInt_t triggerTimeNs = 0;
  UInt_t triggerTimeNsExpected = 0;
  Double_t heading = 0;
  UInt_t l3TrigPatternH = 0;
  Int_t run=0;
  // std::vector<Double_t>* deltaPhiDeg = NULL;

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

  angResTree->Branch("thetaExpected", &thetaExpected);
  angResTree->Branch("phiExpected", &phiExpected);
  angResTree->Branch("triggerTimeNs", &triggerTimeNs);
  angResTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  angResTree->Branch("heading", &heading);
  angResTree->Branch("l3TrigPatternH", &l3TrigPatternH);
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);  
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 5;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);
    // if(header->eventNumber != 60832108) continue;
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){// && header->eventNumber < 60.95e6 && header->eventNumber > 60.85e6){
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
	l3TrigPatternH = header->l3TrigPatternH;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
	globalPhiDeg = usefulEvent->eventNumber;
	
	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();

	cc->correlateEvent(usefulEvent, AnitaPol::kHorizontal);

	TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);

	TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
						       triggeredThetaDeg, l3TrigPatternH);

	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						 zoomThetaDeg, l3TrigPatternH,
						 triggeredPhiDeg, triggeredThetaDeg);
	
	globalPhiDeg = globalPhiDeg < 0 ? globalPhiDeg + 360 : globalPhiDeg;
	globalPhiDeg = globalPhiDeg >= 360 ? globalPhiDeg - 360 : globalPhiDeg;

	triggeredPhiDeg = triggeredPhiDeg < 0 ? triggeredPhiDeg + 360 : triggeredPhiDeg;
	triggeredPhiDeg = triggeredPhiDeg >= 360 ? triggeredPhiDeg - 360 : triggeredPhiDeg;

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
	  delete hTriggeredImageH;
	  delete hGlobalImageH;
	  delete hZoomedImageH;
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
