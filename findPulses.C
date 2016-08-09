#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAdu5Pat.h"
#include "CorrelationSummaryAnita3.h"
#include "CorrelationSummary.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>

const double LONGITUDE_SURF_SEAVEY=167.0564667; ///< Longitude of surface seavey.
const double LATITUDE_SURF_SEAVEY=-77.86185; ///< Latitude of surface seavey.
const double ALTITUDE_SURF_SEAVEY=15.0; ///< Altitude of surface seavey.

const double LONGITUDE_BH=167.06679444; ///< Longitude of borehole antenna.
const double LATITUDE_BH=-77.861936111; ///< Latitude if borehole antenna.
const double ALTITUDE_BH=-33.67; ///< Altitude of borehole antenna.


Long64_t getMinDeltaT(Long64_t x[10], int n);

void findPulses(int firstRun, int lastRun);
  
void findPulses(int firstRun, int lastRun) {

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];

  CalibratedAnitaEvent *event = 0;
  RawAnitaHeader *header =0;
  PrettyAnitaHk *hk = 0;
  Adu5Pat *pat =0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");

  for (int run=firstRun;run<lastRun+1;run++){

  // sprintf(headerName,"/unix/anita3/flight1415/root/run%d/headFile%d.root",run,run);
  // sprintf(hkName,"/unix/anita3/flight1415/root/run%d/prettyHkFile%d.root",run,run);
  // sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsFile%d.root",run,run);
  // sprintf(eventName,"/unix/anita3/calibratedFlight1415/run%d/calEventFile%d.root",run,run);
  sprintf(headerName,"~/UCL/ANITA/flight1415/root/run%d/headFile%d.root",run,run);
  sprintf(hkName,"~/UCL/ANITA/flight1415/root/run%d/prettyHkFile%d.root",run,run);
  sprintf(gpsName,"~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root",run,run);
  sprintf(eventName,"~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root",run,run);

  headChain->Add(headerName);
  gpsChain->Add(gpsName);


  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);

  headChain->BuildIndex("header->eventNumber");
  gpsChain->BuildIndex("pat->realTime");



  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  Int_t cutTimeNs = 60e3;
  char cpol[100];
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  // AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  if (pol == AnitaPol::kVertical){
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    timeOffset = + 92.8;
    //    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }else{ 
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    // sourceLat = - (79 + (27.94097/60)); //WAIS
    // sourceLon = -(112 + (6.76208/60));
    // sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }

  int maxEntry=headChain->GetEntries();

  
  headChain->GetEntry(0);
  UInt_t firstTS = header->realTime;
  headChain->GetEntry(headChain->GetEntries() - 1);
  UInt_t lastTS = header->realTime;


  int tempAnt=-1;

  Double_t deltaTExpected;

  double drealTime;
  double dtriggerTimeNs;
  bool first = true;
  Int_t ndelays = -1;

  UInt_t timedelay[10];
  
  double shortdelay = 3.15e3;

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

  } else {

    timedelay[0] =  50e6;
    timedelay[1] = 250e6;
    timedelay[2] = 450e6;
    timedelay[3] = 650e6;
    timedelay[4] = 850e6;
    ndelays = 5;

    if (firstRun==131 && lastRun==133){
    timedelay[0] =  10e6 +5e6 + 1750;
    timedelay[1] =  60e6 +5e6 + 1750;
    timedelay[2] = 110e6 +5e6 + 1750;
    timedelay[3] = 160e6 +5e6 + 1750;
    timedelay[4] = 210e6 +5e6 + 1750;
    timedelay[5] = 260e6 +5e6 + 1750;
    timedelay[6] = 310e6 +5e6 + 1750;
    timedelay[7] = 360e6 +5e6 + 1750;
    timedelay[8] = 410e6 +5e6 + 1750;
    timedelay[9] = 460e6 +5e6 + 1750;
    ndelays = 10;
    shortdelay = 450.;
    }

  }

  UInt_t timedelay_corr[10];
  for (int i=0;i<ndelays;++i) timedelay_corr[i]= timedelay[i] - i*shortdelay;

  Int_t insertionDelay = 0;//1.1674e3;
  Long64_t deltaTns[10];
  Long64_t deltaTns_corr[10];

  TH1D *h = new TH1D("h", "", 10000, -1e9, 1e9);

  TH1D *hPulse = new TH1D("hPulse", "", 10000, -2e4, 2e4);

  TH1D *hPulse2 = new TH1D("hPulse2", "", 10000, -2e4, 2e4);

  for(Long64_t entry=0;entry<maxEntry;entry++) {
    if(entry%1000==0) std::cout << entry*100./maxEntry << "%        \r" << flush;
    
    headChain->GetEntry(entry);    

    if (pol==AnitaPol::kVertical && header->l3TrigPattern==0) continue;
    if (pol==AnitaPol::kHorizontal && header->l3TrigPatternH==0) continue;
    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);
    
    UsefulAdu5Pat usefulPat(pat);

    UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
    
    UInt_t triggerTimeNs = header->triggerTimeNs;

    h->Fill(int(triggerTimeNsExpected-triggerTimeNs));

    for (int i=0;i<ndelays;++i){
      deltaTns[i]= triggerTimeNsExpected+timedelay[i]+insertionDelay;
      if (deltaTns[i]>1e9) deltaTns[i]-=1e9;
      deltaTns[i] = deltaTns[i] - triggerTimeNs;

      deltaTns_corr[i]= triggerTimeNsExpected+timedelay_corr[i]+insertionDelay;
      if (deltaTns_corr[i]>1e9) deltaTns_corr[i]-=1e9;
      deltaTns_corr[i] = deltaTns_corr[i] - triggerTimeNs;
    }
    
    Long64_t minDeltaT = getMinDeltaT(deltaTns, ndelays);

    hPulse->Fill(minDeltaT);


    Long64_t minDeltaT2 = getMinDeltaT(deltaTns_corr, ndelays);

    hPulse2->Fill(minDeltaT2);

    // if (TMath::Abs(minDeltaT)<60e3) cout << triggerTimeNs << endl;

  }
  


  TCanvas *c1 = new TCanvas("c1", "", 1200, 400);

  c1->SetLogy();

  h->SetTitle(Form("%s pulses; t_{EXPECTED} - t_{MEASURED} (ns);Number of events", cpol));

  h->Draw();


  c1->Print(Form("Pulses_%s_runs_%i_%i.png", cpol, firstRun, lastRun));

  TCanvas *c2 = new TCanvas("c2", "", 1200, 400);

  c2->SetLogy();


  hPulse->SetTitle(Form("%s pulses; t_{EXPECTED}(including delays) - t_{MEASURED} (ns);Number of events", cpol));
  hPulse->Draw();

  c2->Print(Form("Pulses_%s_runs_%i_%i_select.png", cpol, firstRun, lastRun));

  hPulse2->SetTitle(Form("%s pulses; t_{EXPECTED}(including delays) - t_{MEASURED} (ns);Number of events", cpol));
  hPulse2->Draw();

  c2->Print(Form("Pulses_%s_runs_%i_%i_select_corr.png", cpol, firstRun, lastRun));


}



Long64_t getMinDeltaT(Long64_t x[10], int n){

  Long64_t min = 1e20;
  for (int i=0; i<n; i++){
    if ( TMath::Abs(x[i]) < TMath::Abs(min) ){
      min = x[i];
    }
  }
  return min;

}
