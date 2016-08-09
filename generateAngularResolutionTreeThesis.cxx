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

  // // Photogrammetry positions
  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // geom->useKurtAnita3Numbers(1);
  // AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  // for(Int_t surf=0; surf<NUM_SURF; surf++){
  //   for(Int_t chan=0; chan<NUM_CHAN; chan++){
  //     cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
  //   }
  // }

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
  TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);  
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  CrossCorrelator* cc = new CrossCorrelator();
  CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite Notch", 260 - 26, 260 + 26);
  cc->addNotch(notch260);
  CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch", 370 - 26, 370 + 26);
  cc->addNotch(notch370);
  

  
  
  // cc->kZeroChannel16BH = true;
  cc->kUseOffAxisDelay = 0;
  
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
  Double_t zoomPeak2 = 0;
  Double_t zoomPhiDeg2 = 0;
  Double_t zoomThetaDeg2 = 0;
  Double_t deltaThetaDeg2 = 0;
  Double_t deltaPhiDeg2 = 0;
  Double_t zoomPeak3 = 0;
  Double_t zoomPhiDeg3 = 0;
  Double_t zoomThetaDeg3 = 0;
  Double_t deltaThetaDeg3 = 0;
  Double_t deltaPhiDeg3 = 0;
  Double_t zoomPeak4 = 0;
  Double_t zoomPhiDeg4 = 0;
  Double_t zoomThetaDeg4 = 0;
  Double_t deltaThetaDeg4 = 0;
  Double_t deltaPhiDeg4 = 0;
  Double_t zoomPeak5 = 0;
  Double_t zoomPhiDeg5 = 0;
  Double_t zoomThetaDeg5 = 0;
  Double_t deltaThetaDeg5 = 0;
  Double_t deltaPhiDeg5 = 0;
  Double_t zoomPeak6 = 0;
  Double_t zoomPhiDeg6 = 0;
  Double_t zoomThetaDeg6 = 0;
  Double_t deltaThetaDeg6 = 0;
  Double_t deltaPhiDeg6 = 0;
  Double_t zoomPeak7 = 0;
  Double_t zoomPhiDeg7 = 0;
  Double_t zoomThetaDeg7 = 0;
  Double_t deltaThetaDeg7 = 0;
  Double_t deltaPhiDeg7 = 0;
  Double_t zoomPeak8 = 0;
  Double_t zoomPhiDeg8 = 0;
  Double_t zoomThetaDeg8 = 0;
  Double_t deltaThetaDeg8 = 0;
  Double_t deltaPhiDeg8 = 0;

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

  angResTree->Branch("zoomPeak2", &zoomPeak2);
  angResTree->Branch("zoomPhiDeg2", &zoomPhiDeg2);
  angResTree->Branch("zoomThetaDeg2", &zoomThetaDeg2);
  angResTree->Branch("deltaPhiDeg2", &deltaPhiDeg2);
  angResTree->Branch("deltaThetaDeg2", &deltaThetaDeg2);

  angResTree->Branch("zoomPeak3", &zoomPeak3);
  angResTree->Branch("zoomPhiDeg3", &zoomPhiDeg3);
  angResTree->Branch("zoomThetaDeg3", &zoomThetaDeg3);
  angResTree->Branch("deltaPhiDeg3", &deltaPhiDeg3);
  angResTree->Branch("deltaThetaDeg3", &deltaThetaDeg3);

  angResTree->Branch("zoomPeak4", &zoomPeak4);
  angResTree->Branch("zoomPhiDeg4", &zoomPhiDeg4);
  angResTree->Branch("zoomThetaDeg4", &zoomThetaDeg4);
  angResTree->Branch("deltaPhiDeg4", &deltaPhiDeg4);
  angResTree->Branch("deltaThetaDeg4", &deltaThetaDeg4);

  angResTree->Branch("zoomPeak5", &zoomPeak5);
  angResTree->Branch("zoomPhiDeg5", &zoomPhiDeg5);
  angResTree->Branch("zoomThetaDeg5", &zoomThetaDeg5);
  angResTree->Branch("deltaPhiDeg5", &deltaPhiDeg5);
  angResTree->Branch("deltaThetaDeg5", &deltaThetaDeg5);

  angResTree->Branch("zoomPeak6", &zoomPeak6);
  angResTree->Branch("zoomPhiDeg6", &zoomPhiDeg6);
  angResTree->Branch("zoomThetaDeg6", &zoomThetaDeg6);
  angResTree->Branch("deltaPhiDeg6", &deltaPhiDeg6);
  angResTree->Branch("deltaThetaDeg6", &deltaThetaDeg6);

  angResTree->Branch("zoomPeak7", &zoomPeak7);
  angResTree->Branch("zoomPhiDeg7", &zoomPhiDeg7);
  angResTree->Branch("zoomThetaDeg7", &zoomThetaDeg7);
  angResTree->Branch("deltaPhiDeg7", &deltaPhiDeg7);
  angResTree->Branch("deltaThetaDeg7", &deltaThetaDeg7);

  angResTree->Branch("zoomPeak8", &zoomPeak8);
  angResTree->Branch("zoomPhiDeg8", &zoomPhiDeg8);
  angResTree->Branch("zoomThetaDeg8", &zoomThetaDeg8);
  angResTree->Branch("deltaPhiDeg8", &deltaPhiDeg8);
  angResTree->Branch("deltaThetaDeg8", &deltaThetaDeg8);
  
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
  Double_t maxDPhiDeg2;
  angResTree->Branch("maxDPhiDeg2", &maxDPhiDeg2);
  Double_t maxDPhiDeg3;
  angResTree->Branch("maxDPhiDeg3", &maxDPhiDeg3);
  Double_t maxDPhiDeg4;
  angResTree->Branch("maxDPhiDeg4", &maxDPhiDeg4);
  Double_t maxDPhiDeg5;
  angResTree->Branch("maxDPhiDeg5", &maxDPhiDeg5);
  Double_t maxDPhiDeg6;
  angResTree->Branch("maxDPhiDeg6", &maxDPhiDeg6);
  Double_t maxDPhiDeg7;
  angResTree->Branch("maxDPhiDeg7", &maxDPhiDeg7);
  Double_t maxDPhiDeg8;
  angResTree->Branch("maxDPhiDeg8", &maxDPhiDeg8);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 115;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    if(60841774==header->eventNumber){
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
	  thetaExpected*=TMath::RadToDeg();
	
	  phiExpected = phiExpected < 0 ? phiExpected + 360 : phiExpected;
	  phiExpected = phiExpected >= 360 ? phiExpected - 360 : phiExpected;

	  const int ant1 = 14;
	  // const int ant2 = 31;

	  int p1 = ant1 % NUM_PHI;
	  int r1 = ant1 / NUM_PHI;
    
	  TString label1 = TString::Format("%d", p1+1);
	  if(r1==0){
	    label1 += "TH";
	  }
	  else if(r1==1){
	    label1 += "MH";
	  }
	  else{
	    label1 += "BH";
	  }

	  const int numAnt2s = 7;
	  int theAnt2s[numAnt2s] = {-1, 15, 31, 47, 0, 16, 32};

	  for(int ant2Ind = 0; ant2Ind < numAnt2s; ant2Ind++){
	    const int ant2 = theAnt2s[ant2Ind];
	    const int comboInd = ant2 >= 0 ? cc->comboIndices[ant1][ant2] : -1;
	    
	    int p2 = ant2 % NUM_PHI;
	    int r2 = ant2 / NUM_PHI;
    
	    TString label2 = TString::Format("%d", p2+1);
	    if(r2==0){
	      label2 += "TH";
	    }
	    else if(r2==1){
	      label2 += "MH";
	    }
	    else{
	      label2 += "BH";
	    }

	    
	    cc->kOnlyThisCombo = comboInd;
	    std::cout << comboInd << std::endl;
	    cc->eventNumber[0] = 0; // force reevaltulation of the cross-correlations
	    // cc->lastEventUpsampleCorrelated[0] = 0;	    

	    cc->correlateEvent(usefulEvent, pol);

	    if(ant2Ind==0){
	      TGraph* gr1 = cc->grs[AnitaPol::kHorizontal][ant1];
	      TGraph* gr1I = cc->grsResampled[AnitaPol::kHorizontal][ant1];
	      for(int i=0; i < gr1I->GetN(); i++){
		gr1I->GetY()[i]*=cc->interpRMS[AnitaPol::kHorizontal][ant1];
	      }
	      gr1->Write();
	      gr1I->Write();
	      // delete gr1;
	      // delete gr1I;
	    }

	    if(comboInd >= 0){
	      TGraph* gr2 = cc->grs[AnitaPol::kHorizontal][ant2];
	      TGraph* gr2I = cc->grsResampled[AnitaPol::kHorizontal][ant2];
	      for(int i=0; i < gr2I->GetN(); i++){
		gr2I->GetY()[i]*=cc->interpRMS[AnitaPol::kHorizontal][ant2];
	      }
	      gr2->Write();
	      gr2I->Write();
	      // delete gr2;
	      // delete gr2I;
	    
	      TGraph* grCor = cc->getCrossCorrelationGraph(AnitaPol::kHorizontal, ant1, ant2);
	      grCor->Write();
	      delete grCor;
	    }	    
	  	
	    TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);

	    if(comboInd >= 0){
	      TString title = TString::Format("Event number %u - Cross-correlation %s & %s",
					      eventNumber,
					      label1.Data(),
					      label2.Data());
	      hGlobalImageH->SetTitle(title);
	      hGlobalImageH->SetName(TString::Format("hGlobalImageH_%u_%d", eventNumber, ant2));	      
	    }
	    else{
	      hGlobalImageH->Write();
	      delete hGlobalImageH;
	      hGlobalImageH = NULL;
	    }
	    
	    TH2D* hZoomedImageH = cc->makeZoomedImage(AnitaPol::kHorizontal,
						      zoomPeak, zoomPhiDeg, zoomThetaDeg,
						      globalPhiDeg, globalThetaDeg);

	    if(comboInd >= 0){

	      hZoomedImageH->SetName(TString::Format("hZoomedImageH_%u_%d", eventNumber, ant2));
	      TGraph* grCor2 = cc->getUpsampledCrossCorrelationGraph(AnitaPol::kHorizontal, ant1, ant2);
	      grCor2->Write();
	      delete grCor2;
	    }  
	    else{
	      hZoomedImageH->Write();
	      delete hZoomedImageH;
	      hZoomedImageH = NULL;
	    }

	    // TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
	    // 						    triggeredThetaDeg, l3TrigPatternH);
	
	    if(numSaved < maxToSave){
	      numSaved++;
	    }
	    else{
	      delete hGlobalImageH;
	      delete hZoomedImageH;
	    }

	  }

	  

	  angResTree->Fill();

	  delete usefulEvent;
	}
      }
    }
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
