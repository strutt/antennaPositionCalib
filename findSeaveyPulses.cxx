// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to look at runs and figure out how far we were from LDB.
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"

#include "ProgressBar.h"

const Double_t sourceLat = - (77 + (51.23017/60));
const Double_t sourceLon = +(167 + (12.16908/60));
const Double_t sourceAlt = 0;

int main(int argc, char *argv[])
{
  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 0;
  }

  // 145
  // 161
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");

  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);

  for(int run=firstRun; run<lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
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

  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2000;  
  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");    
  TString title1 = TString::Format("Distance from LDB runs %d - %d; run; Distance (km)", firstRun, lastRun);
  TH2D* hDistanceFromLdb = new TH2D("hDistanceFromLdb", title1,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDistBins, 0, maxDistKm);


  const Double_t minDiff = 24.994e6;
  const Double_t maxDiff = 25.003e6;
  // const Double_t minDiff = 24e6; //24.985e6;
  // const Double_t maxDiff = 76e6; //25.005e6;  
  const Int_t numDiffBins = 1024;
  TString title2 = TString::Format("Event TriggerTimeNs difference from expected LDB arrival time runs %d - %d; run; #Delta t_{received - expected} (ns)", firstRun, lastRun);
  TH2D* hDeltaTFromLdb = new TH2D("hTriggerTimeNs", title2,
				     lastRun+1-firstRun, firstRun, lastRun+1, numDiffBins, minDiff, maxDiff);
  
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){
    
      UsefulAdu5Pat usefulPat(pat);
      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
      hDistanceFromLdb->Fill(header->run, distKm);

      Int_t a = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected) - 0.0005e6; // + 0.003e8;
      a = a%Int_t(1999969e2);
      
      if(a > 67.5e6){
      	a -= 50.0012e6;
      }
      else if(a > 27.5e6){
      	a -= 25e6;
      }
      // std::cout << a << "\t" << a%Int_t(2e8) << std::endl;
      hDeltaTFromLdb->Fill(header->run, a);
    
      // std::cout << distKm << std::endl;

    }
    p++;
  }


  outFile->Write();
  outFile->Close();

  return 0;

}
