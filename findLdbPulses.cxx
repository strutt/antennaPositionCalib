// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to look at runs and figure out how far we were from WAIS divide.
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

#include <ProgressBar.h>

int main(int argc, char *argv[])
{
  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 0;
  }

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");

  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);
  // const Int_t cutTimeNs = Int_t(10e3);
  // const Int_t cutTimeNs = Int_t(60e3);  
  // const Int_t cutTimeNs = Int_t(60e3);  
  const Int_t cutTimeNs = Int_t(3000);
  
  for(int run=firstRun; run<lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);    
    gpsChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");  
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //5000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");    


  const Double_t sourceLat = - (77 + (51.23017/60)); // OLD NUMBERS
  // const Double_t sourceLat = - (77 + 51./60 + 44.16/3600); // From Google Earth 77°51'44.16"S
  const Double_t sourceLon = +(167 + (12.16908/60));// OLD NUMBERS
  // const Double_t sourceLon = +(167 + 2./60 + 28.83/3600); // From Google Earth 167° 2'28.83"E
  const Double_t sourceAlt = 0;
  

  headChain->GetEntry(0);
  const int firstEventNumber = header->eventNumber;
  headChain->GetEntry(maxEntry-1);
  const int lastEventNumber = header->eventNumber;

  // const Double_t maxDeltaTriggerTimeNs = cutTimeNs;
  const Double_t maxDeltaTriggerTimeNs = cutTimeNs;  
  const Double_t minDiff = -maxDeltaTriggerTimeNs;
  const Double_t maxDiff =  maxDeltaTriggerTimeNs;
  // const Double_t minDiff = -1e6;
  // const Double_t maxDiff = 1e6;


  UInt_t timedelay[10];
  Int_t shortdelay = 3150;
  int ndelays = -1;
  Int_t deltaTns_corr[5];

  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2000;
  const int nEventBins = 1024;
  TString title1 = TString::Format("Distance from LDB runs %d - %d; Event Number; Distance (km)", firstRun, lastRun);
  TH2D* hDistanceFromWais = new TH2D("hDistanceFromWais", title1,
				     nEventBins, firstEventNumber, lastEventNumber+1,
				     numDistBins, 0, maxDistKm);

  const Int_t numDiffBins = 128;

  TString title2 = TString::Format("#Deltat_{trigger-expected} difference from expected LDB runs %d - %d; Event Number; #Delta t_{trigger - expected} (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromWais = new TH2D("hTriggerTimeNs", title2,
				     nEventBins, firstEventNumber, lastEventNumber+1, numDiffBins, minDiff, maxDiff);



  TString title3 = TString::Format("TriggerTimeNs for RF triggered events runs %d - %d; Event Number; triggerTimeNs (ns)", firstRun, lastRun);  
  TH2D* hDeltaTFromWais2 = new TH2D("hTriggerTimeNs2", title3,
				    nEventBins, firstEventNumber, lastEventNumber+1, numDiffBins, 0, 1e9);
  
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    // gpsChain->GetEntryWithIndex(header->realTime);
    gpsChain->GetEntry(entry);

    if (header->run<=153){
      timedelay[0] = Int_t(25e6);
      timedelay[1] = Int_t(225e6);
      timedelay[2] = Int_t(425e6);
      timedelay[3] = Int_t(625e6);
      timedelay[4] = Int_t(825e6); // in ms
    }
    else{//      if (firstRun>153){
      timedelay[0] = Int_t(50e6);
      timedelay[1] = Int_t(250e6);
      timedelay[2] = Int_t(450e6);
      timedelay[3] = Int_t(650e6);
      timedelay[4] = Int_t(850e6);
    }
    ndelays = 5;
  
    Long64_t timedelay_corr[10];
    for (int i=0;i<ndelays;++i){
      timedelay_corr[i]= timedelay[i] - i*shortdelay;
    }
  
    if((header->trigType & 1)==1){
    
      UsefulAdu5Pat usefulPat(pat);
      // UInt_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      UInt_t triggerTimeNs = header->triggerTimeNs;
      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt); 
      Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
      hDistanceFromWais->Fill(header->eventNumber, distKm);

      Long64_t minDeltaT = Long64_t(1e11);
      for (int i=0;i<ndelays;++i){
	deltaTns_corr[i]= triggerTimeNsExpected+timedelay_corr[i];
	if(deltaTns_corr[i]>Int_t(1e9)){
	  deltaTns_corr[i]-=Int_t(1e9);
	}
	deltaTns_corr[i] = deltaTns_corr[i] - triggerTimeNs;
	if(TMath::Abs(deltaTns_corr[i]) < TMath::Abs(minDeltaT) ){
	  minDeltaT = deltaTns_corr[i];
	}
      }
      if(TMath::Abs(minDeltaT) < cutTimeNs){
	// if(header->run==152){
	//   std::cout << header->eventNumber << std::endl;
	// }
	
	// Int_t deltaTriggerTimeNs = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected);
	hDeltaTFromWais->Fill(header->eventNumber, minDeltaT);
    
      // std::cout << distKm << std::endl;
      }
      hDeltaTFromWais2->Fill(header->eventNumber, triggerTimeNs);      
    }
    p++;
  }


  outFile->Write();
  outFile->Close();

  return 0;

}
