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

  TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, AnitaPol::kVertical);    
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  CrossCorrelator* cc = new CrossCorrelator();
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

  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  lindaFileNameReference->Write();
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //1000; //10000;
  Long64_t startEntry = 0; //
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    if((header->trigType & 1)==1){// && numPhiTrigs==2){// && trigIWant > 0){
      gpsChain->GetEntry(entry);
      
      UsefulAdu5Pat usefulPat(pat);
      Double_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      Double_t triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){

	calEventChain->GetEntry(entry);

	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	Double_t phiExpected;
	Double_t thetaExpected;
  	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();

	cc->correlateEvent(usefulEvent, pol);

	// TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);

	UInt_t l3TrigPatternH = header->l3TrigPatternH;
	Double_t zoomPeak, zoomPhiDeg, zoomThetaDeg;
	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						 zoomThetaDeg, l3TrigPatternH,
						 phiExpected, thetaExpected);

	
	delete usefulEvent;
      }	
    }    
    p++;
  }
  
  outFile->Write();
  outFile->Close();
  return 0;
}


