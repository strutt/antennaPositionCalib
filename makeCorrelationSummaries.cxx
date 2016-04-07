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
 
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //1000; //1000; //10000;
  Long64_t startEntry = 0; //
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  std::map<UInt_t, Int_t> count;
  std::map<UInt_t, TH2D*> hAveSummary;

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    // if(header->eventNumber != 60832108) continue;

    UInt_t l3TrigPatternH = 0;
    l3TrigPatternH = header->l3TrigPatternH;

    Int_t numPhiTrigs = 0;
    for(int phi=0; phi<NUM_PHI; phi++){
      Int_t state = RootTools::getBit(phi, l3TrigPatternH);
      if(state > 0){
	numPhiTrigs++;
      }
    }

    if((header->trigType & 1)==1 && numPhiTrigs==2){// && trigIWant > 0){
      gpsChain->GetEntry(entry);
      
      UsefulAdu5Pat usefulPat(pat);
      Double_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      Double_t triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){

	calEventChain->GetEntry(entry);

	UInt_t l3TrigPattern = 0;
	l3TrigPattern = header->l3TrigPattern;

	
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	Double_t phiExpected;
	Double_t thetaExpected;	
	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);


	phiExpected*=TMath::RadToDeg();
	thetaExpected*=TMath::RadToDeg();

	cc->correlateEvent(usefulEvent, pol);
    
	Double_t triggeredPeak = 0;
	Double_t triggeredPhiDeg = 0;
	Double_t triggeredThetaDeg = 0;
	Double_t zoomPeak = 0;
	Double_t zoomPhiDeg = 0;
	Double_t zoomThetaDeg = 0;

	TH2D* hTriggeredImageH2 = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
							 triggeredThetaDeg, l3TrigPatternH);

	TH2D* hZoomedImageH2 = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						   zoomThetaDeg, l3TrigPatternH,
						   triggeredPhiDeg, triggeredThetaDeg);

	// TH2D* hSummary = cc->makeCorrelationSummaryHistogram(pol, l3TrigPatternH, zoomPhiDeg, zoomThetaDeg);
	TH2D* hSummary2 = cc->makeCorrelationSummaryHistogram(pol, l3TrigPatternH, zoomPhiDeg, zoomThetaDeg);
	hSummary2->SetName("hTemp");


	TH2D* hSummary = cc->makeCorrelationSummaryHistogram(pol, l3TrigPatternH, phiExpected, thetaExpected);

	// Now lets try looking at the correlation values at the peak vs. the expected peak.
	for(int biny=1; biny<=hSummary->GetNbinsY(); biny++){
	  for(int binx=1; binx<=hSummary->GetNbinsX(); binx++){
	    Double_t val1 = hSummary->GetBinContent(binx, biny);
	    Double_t val2 = hSummary2->GetBinContent(binx, biny);
	    hSummary->SetBinContent(binx, biny, val1-val2);
	  }
	}

	delete hSummary2;
	
	if(count.count(l3TrigPatternH)==0){
	  hAveSummary[l3TrigPatternH] = hSummary;
	  hSummary->SetName(TString::Format("h%u", l3TrigPatternH));
	  hSummary->SetTitle(TString::Format("Average correlation at peak value for l3TrigPattern = %u", l3TrigPatternH));
	  count[l3TrigPatternH] = 1;
	}
	else{
	  TH2D* h = hAveSummary[l3TrigPatternH];
	  for(int biny=1; biny<=h->GetNbinsY(); biny++){
	    for(int binx=1; binx<=h->GetNbinsX(); binx++){
	      Double_t val1 = hSummary->GetBinContent(binx, biny);
	      Double_t val2 = h->GetBinContent(binx, biny);	      
	      h->SetBinContent(binx, biny, val1+val2);
	    }
	  }
	  count[l3TrigPatternH]++;
	  delete hSummary;
	}
    
	delete hTriggeredImageH2;
	delete hZoomedImageH2;


	delete usefulEvent;

      }
    }
    p++;
  }

  std::map<UInt_t, TH2D*>::iterator i;
  for(i = hAveSummary.begin(); i != hAveSummary.end(); ++i){
    TH2D* h = (*i).second;
    Int_t theCount = count[(*i).first];
    for(int biny=1; biny<=h->GetNbinsY(); biny++){
      for(int binx=1; binx<=h->GetNbinsX(); binx++){
	Double_t sum = h->GetBinContent(binx, biny);
	h->SetBinContent(binx, biny, sum/theCount);
      }
    }
    TNamed* n = new TNamed(TString::Format("n%u", (*i).first),
			   TString::Format("%d", theCount));
    n->Write();
    delete n;
  }
  
  outFile->Write();
  outFile->Close();
  return 0;
}


