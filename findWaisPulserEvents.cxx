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


  headChain->GetEntry(0);
  const int firstEventNumber = header->eventNumber;
  headChain->GetEntry(maxEntry-1);
  const int lastEventNumber = header->eventNumber;

  const Double_t maxDeltaTriggerTimeNs = 1200;
  const Double_t minDiff = -2*maxDeltaTriggerTimeNs;
  const Double_t maxDiff =  2*maxDeltaTriggerTimeNs;
  // const Double_t minDiff = -1e6;
  // const Double_t maxDiff = 1e6;

  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2000;
  const int nEventBins = 1024;
  TString title1 = TString::Format("Distance from WAIS divide runs %d - %d; Event Number; Distance (km)", firstRun, lastRun);
  TH2D* hDistanceFromWais = new TH2D("hDistanceFromWais", title1,
				     nEventBins, firstEventNumber, lastEventNumber+1,
				     numDistBins, 0, maxDistKm);

  const Int_t numDiffBins = 128;

  TString title2 = TString::Format("Event TriggerTimeNs difference from expected WAIS arrival time runs %d - %d; Event Number; #Delta t_{trigger - expected} (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromWais = new TH2D("hTriggerTimeNs", title2,
				     nEventBins, firstEventNumber, lastEventNumber+1, numDiffBins, minDiff, maxDiff);


  TString title3 = TString::Format("TriggerTimeNs for RF triggered events runs %d - %d; Event Number; triggerTimeNs (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromWais2 = new TH2D("hTriggerTimeNs2", title3,
				     nEventBins, firstEventNumber, lastEventNumber+1, numDiffBins, 1e6, 2e6);
  
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    // gpsChain->GetEntryWithIndex(header->realTime);
    gpsChain->GetEntry(entry);    
    if((header->trigType & 1)==1){
    
      UsefulAdu5Pat usefulPat(pat);
      UInt_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
      hDistanceFromWais->Fill(header->eventNumber, distKm);

      Int_t deltaTriggerTimeNs = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected);
      hDeltaTFromWais->Fill(header->eventNumber, deltaTriggerTimeNs);
    
      // std::cout << distKm << std::endl;

    }
    hDeltaTFromWais2->Fill(header->eventNumber, header->triggerTimeNs);
    
    p++;
  }


  outFile->Write();
  outFile->Close();

  return 0;

}
