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

  // Validate Linda's Numbers
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

  TString lindaFileName = "newLindaNumbers_4steps_VPOL_10kVSeavey_2015_12_17_time_17_41_28";

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


  CrossCorrelator* cc = new CrossCorrelator();  
  
  // Double_t sourceLat = - (77 + (51.23017/60));
  // Double_t sourceLon = +(167 + (12.16908/60));
  // Double_t sourceAlt = 0;
  const Double_t sourceLat = - (77 + (51.23017/60));
  const Double_t sourceLon = +(167 + (12.16908/60));
  const Double_t sourceAlt = 0;

  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;
  const Int_t cutTimeNs = 60e3;
  // const Double_t minDeltaTriggerTimeNs = 24.985e6;
  // const Double_t maxDeltaTriggerTimeNs = 25.005e6;
  // const Double_t minDeltaTriggerTimeNs = 24.994e6;
  // const Double_t maxDeltaTriggerTimeNs = 25.003e6;

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // geom->useKurtAnita3Numbers(1);
  // AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  // for(int surf=0; surf<NUM_SURF; surf++){
  //   for(int chan=0; chan<NUM_CHAN; chan++){
  //     cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
  //   }
  // }

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("/unix/anita3/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("/unix/anita3/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("/unix/anita3/flight1415/root/run%d/calEventFile%d.root", run, run);
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
  UInt_t l3TrigPattern = 0;
  Int_t run = 0;  
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
  angResTree->Branch("l3TrigPattern", &l3TrigPattern);
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);  

  UInt_t timedelay[10];
  Int_t shortdelay = 3150;
  int ndelays = -1;
  Int_t deltaTns_corr[5];
  if (pol==AnitaPol::kVertical){
    timedelay[0] = 25e6;
    timedelay[1] = 225e6;
    timedelay[2] = 425e6;
    timedelay[3] = 625e6;
    timedelay[4] = 825e6; // in ms
    if (firstRun>153){
      timedelay[0] =  50e6;
      timedelay[1] = 250e6;
      timedelay[2] = 450e6;
      timedelay[3] = 650e6;
      timedelay[4] = 850e6;
    }
    ndelays = 5;

  } 

  Long64_t timedelay_corr[10];
  for (int i=0;i<ndelays;++i) timedelay_corr[i]= timedelay[i] - i*shortdelay;


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
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){// && header->eventNumber < 60.95e6 && header->eventNumber > 60.85e6){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      triggerTimeNs = header->triggerTimeNs;

      Long64_t minDeltaT = 1e11;
      for (int i=0;i<ndelays;++i){
	deltaTns_corr[i]= triggerTimeNsExpected+timedelay_corr[i];
	if (deltaTns_corr[i]>1e9) deltaTns_corr[i]-=1e9;
	deltaTns_corr[i] = deltaTns_corr[i] - triggerTimeNs;
	if ( TMath::Abs(deltaTns_corr[i]) < TMath::Abs(minDeltaT) ){
	  minDeltaT = deltaTns_corr[i];
	}
      }
      
     if(TMath::Abs(minDeltaT) < cutTimeNs){

      // Int_t a = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected) - 0.0005e6; // + 0.003e8;
      // a = a%Int_t(1999969e2);
      // if(a > 67.5e6){
      // 	a -= 50.0012e6;
      // }
      // else if(a > 27.5e6){
      // 	a -= 25e6;
      // }
      // if(a > minDeltaTriggerTimeNs && a < maxDeltaTriggerTimeNs){
      // if(true){
	eventNumber = header->eventNumber;
	run = header->run;
	// std::cout << "eventNumber " << header->eventNumber << std::endl;

	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	l3TrigPattern = header->l3TrigPattern;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
	globalPhiDeg = usefulEvent->eventNumber;
	
	// usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	// usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=TMath::RadToDeg();

	cc->correlateEvent(usefulEvent, pol);

	TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);

	TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
						       triggeredThetaDeg, l3TrigPattern);

	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						 zoomThetaDeg, l3TrigPattern,
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
