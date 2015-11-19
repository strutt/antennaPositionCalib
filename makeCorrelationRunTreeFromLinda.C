#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "CorrelationSummaryAnita3.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

void makeCorrelationRunTree(int run, int numEnts=0, char *outDir=0);

Long64_t getMinDeltaT(Long64_t x[5]);

void makeCorrelationRunTree(int run, int numEnts, char *outDir) {
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

   AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

   Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
   Int_t cutTimeNs = 60e3;
   char cpol[100];


   if (pol == AnitaPol::kVertical){
     sourceLat = - (77 + (51.23017/60));
     sourceLon = +(167 + (12.16908/60));
     sourceAlt = 0;
     //     timeOffset = + 92.8;
    sprintf(cpol, "VPOL");
  }else{ 
    sourceLat = - (79 + (27.94097/60));
    sourceLon = -(112 + (6.76208/60));
    sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }


  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];

  CalibratedAnitaEvent *event = NULL;
  RawAnitaHeader *header =NULL;
  PrettyAnitaHk *hk = NULL;
  Adu5Pat *pat = NULL;
  
  TChain *eventChain = new TChain("eventTree");
  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  TChain *prettyHkChain = new TChain("prettyHkTree");

//   for (int run=firstRun;run<lastRun+1;run++){
  //sprintf(eventName,"/unix/anita3/flight1415/root/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/unix/anita3/flight1415/root/run%d/headFile%d.root",run,run);
  sprintf(hkName,"/unix/anita3/flight1415/root/run%d/prettyHkFile%d.root",run,run);
  sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsFile%d.root",run,run);
  sprintf(eventName,"/unix/anita3/flight1415/root/run%d/calEventFile%d.root",run,run);

  eventChain->Add(eventName);
  headChain->Add(headerName);
  prettyHkChain->Add(hkName);
  gpsChain->Add(gpsName);

//   }
  eventChain->SetBranchAddress("event",&event);
  headChain->SetBranchAddress("header",&header);
  prettyHkChain->SetBranchAddress("hk",&hk);
  gpsChain->SetBranchAddress("pat",&pat);

//   eventChain->BuildIndex("event->eventNumber");
  prettyHkChain->BuildIndex("hk->realTime");
  gpsChain->BuildIndex("pat->realTime");


  //Make output files
  CorrelationSummaryAnita3 *theCor=0;
  TFile *fpOut;
  sprintf(outName,"%s/corRun_NEW3_%s_seavey_%d.root",outDir,cpol, run);
  
  cout << outName << endl;
  fpOut= new TFile(outName,"RECREATE");
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummaryAnita3",&theCor);

  Long64_t maxEntry=headChain->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/1000;
  if(starEvery==0) starEvery=1;
  
  Double_t thetaWave,phiWave;

  UInt_t triggerTimeNsExpected;
  UInt_t triggerTimeNs;
  int ant;

  cout << maxEntry << endl;

  Double_t deltaT= 1. / (2.6);

  double c3poNumStart = 0;
  int count = 0;
  int corrected = 0;

  double max =100;

  double dc3poNum, dtrigTime;
  
  UInt_t timedelay[5] = {25e6, 225e6, 425e6, 625e6, 825e6}; // in ms
  
  if (run>153){
    timedelay[0] =  50e6;
    timedelay[1] = 250e6;
    timedelay[2] = 450e6;
    timedelay[3] = 650e6;
    timedelay[4] = 850e6;
  }
 
  UInt_t timedelay_corr[5];
  for (int i=0;i<5;++i) timedelay_corr[i]= timedelay[i] - i*3.15e3;


  Int_t insertionDelay = 0;//1.1674e3;
  Long64_t deltaTns[5];
  Long64_t deltaTns_corr[5];

  for(Long64_t entry=0;entry<maxEntry;entry++) {
//  for(Long64_t entry=0;entry<10000;entry++) {

     if (entry%5000==0) cout << entry*100./maxEntry << " %               \r" << flush;
    
    headChain->GetEntry(entry);
    
    if (header->l3TrigPattern==0) continue;

    eventChain->GetEntry(entry);    

    Long64_t prettyHkEntry = prettyHkChain->GetEntryNumberWithIndex(header->realTime);
    if(prettyHkEntry < 0 ) continue;
    prettyHkChain->GetEntry(prettyHkEntry);
    //     prettyHkChain->GetEntry(entry);

    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);
        

     UsefulAdu5Pat usefulPat(pat);

     triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

     triggerTimeNs = Double_t(header->triggerTimeNs);
     
    for (int i=0;i<5;++i){
      deltaTns[i]= triggerTimeNsExpected+timedelay[i]+insertionDelay;
      if (deltaTns[i]>1e9) deltaTns[i]-=1e9;
      deltaTns[i] = deltaTns[i] - triggerTimeNs;

      deltaTns_corr[i]= triggerTimeNsExpected+timedelay_corr[i]+insertionDelay;
      if (deltaTns_corr[i]>1e9) deltaTns_corr[i]-=1e9;
      deltaTns_corr[i] = deltaTns_corr[i] - triggerTimeNs;
    }
    
    Long64_t minDeltaT = getMinDeltaT(deltaTns_corr);


     if(TMath::Abs(minDeltaT) > cutTimeNs) continue;

     PrettyAnitaEvent realEvent(event,WaveCalType::kDefault,hk);

     usefulPat.getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,thetaWave,phiWave);
     ant=fGeomTool->getTopAntNearestPhiWave(phiWave, pol);

//      int maxAnt = realEvent.getMaxAntenna(pol);
//      while (maxAnt>15) maxAnt-=16;
//      cout << ant << " " << maxAnt << endl;

     theCor =realEvent.getCorrelationSummaryAnita3(ant,pol,deltaT);
     corTree->Fill();     
     delete theCor;
     
  }
  std::cerr << "\n";
  corTree->AutoSave();
  fpOut->Close();
}
     
Long64_t getMinDeltaT(Long64_t x[5]){

  Long64_t min = 1e20;

  for (int i=0; i<5; i++){
    if ( TMath::Abs(x[i]) < TMath::Abs(min) ){
      min = x[i];
    }
  }
  return min;

}
