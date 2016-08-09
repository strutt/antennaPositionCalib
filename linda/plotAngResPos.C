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

void plotAngResPos(string arSample, string pcSample);

void plotAngResPos(string arSample, string pcSample) {
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

   // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
   AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
   string folderName;
   int firstRun=-1;
   int lastRun=-1;
   string baseName = "/unix/anita3/strutt/antennaPositionCalibPlots/newFittingOrder/";

  if (pcSample=="WAIS"){
    //    folderName = "newLindaNumbers_4steps_WAISHPOL_2015_11_19_time_15_06_17";
    folderName = "../ldbPulses/newLindaNumbers_4steps_WAISHPOL_2015_11_27_time_12_13_09";
  } else if (pcSample=="LDBHPOL") {
    folderName = "newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45";
  } else if (pcSample=="LDBVPOL") {
    folderName = "newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_19_time_11_49_04";
  }

  // choose the sample for the angular resolution plots
  if (arSample=="WAIS"){
    firstRun = 332;
    lastRun = 354;
  } else if (arSample=="LDBHPOL") {
    firstRun = 151;
    lastRun = 153;
  } else if (arSample=="LDBVPOL") {
    firstRun = 145;
    lastRun = 161;
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
  
  TChain *gpsChain = new TChain("adu5PatTree");
  //  TChain *eventChain = new TChain("eventTree");
  TChain *headChain = new TChain("headTree");
  //  TChain *prettyHkChain = new TChain("prettyHkTree");

  TChain *resChain = new TChain("angResTree");

  for (int irun=firstRun; irun<lastRun+1; irun++){
    cout << Form("%s%s/*%d*root", baseName.c_str(), folderName.c_str(), irun)<< endl;
    resChain->Add(Form("%s%s/*%d*17-41*root", baseName.c_str(), folderName.c_str()), irun);

    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/headFile%d.root",irun,irun);
    sprintf(hkName,"/unix/anita3/flight1415/root/run%d/prettyHkFile%d.root",irun,irun);
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsFile%d.root",irun,irun);
    sprintf(eventName,"/unix/anita3/flight1415/root/run%d/calEventFile%d.root",irun,irun);

    //    eventChain->Add(eventName);
    headChain->Add(headerName);
    //    prettyHkChain->Add(hkName);
    gpsChain->Add(gpsName);

  }
  //  eventChain->SetBranchAddress("event",&event);
  headChain->SetBranchAddress("header",&header);
  //  prettyHkChain->SetBranchAddress("hk",&hk);
  gpsChain->SetBranchAddress("pat",&pat);

  headChain->BuildIndex("header->eventNumber");
  //  eventChain->BuildIndex("event->eventNumber");
  //  prettyHkChain->BuildIndex("hk->realTime");
  gpsChain->BuildIndex("pat->realTime");

  UInt_t eventNumber;
  Double_t deltaPhiDeg;
  Double_t deltaThetaDeg;
  resChain->SetBranchAddress("eventNumber", &eventNumber);
  resChain->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);
  resChain->SetBranchAddress("deltaThetaDeg", &deltaThetaDeg);

  //Make output files

  Long64_t maxEntry=resChain->GetEntries(); 

  Int_t starEvery=maxEntry/1000;
  if(starEvery==0) starEvery=1;
  
  Double_t thetaWave,phiWave;

  std::cout << "All entries " << maxEntry << std::endl;

  Double_t deltaT= 1. / (2.6*40.);


  headChain->GetEntry(0);
  double firstTS = header->realTime;
  headChain->GetEntry(headChain->GetEntries()-1);
  double lastTS = header->realTime;


  TH2D *hPhi_time = new TH2D("hPhi_time", "", 100, firstTS, lastTS, 100, -2, 2);
  TH2D *hTheta_time = new TH2D("hTheta_time", "", 200, firstTS, lastTS, 100, -2, 2);

  cout << firstTS-1.41894e+09 << " " << lastTS-1.41894e+09 << endl;
  headChain->GetEntry(0);

  
  for(Long64_t entry=0;entry<maxEntry;entry++) {
//  for(Long64_t entry=0;entry<10000;entry++) {

//     if (entry%5000==0) cout << entry*100./maxEntry << " %               \r" << endl;

    resChain->GetEntry(entry);
    Long64_t headEntry = headChain->GetEntryNumberWithIndex(eventNumber);
    // cout << eventNumber << " " << headEntry << endl;

    // break;
    if (headEntry<0) continue;
    headChain->GetEntry(headEntry);


    //    Long64_t eventEntry = eventChain->GetEntryNumberWithIndex(eventNumber);
    //    eventChain->GetEntry(eventEntry);    
    //    Long64_t prettyHkEntry = prettyHkChain->GetEntryNumberWithIndex(header->realTime);
    //    if(prettyHkEntry < 0 ) continue;
    //    prettyHkChain->GetEntry(prettyHkEntry);
    //     prettyHkChain->GetEntry(entry);

    // Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    // if(gpsEntry < 0 ) continue;
    // gpsChain->GetEntry(gpsEntry);
        

    //  UsefulAdu5Pat usefulPat(pat);

     //     PrettyAnitaEvent realEvent(event,WaveCalType::kDefault,hk);


    // cout << header->realTime << " " << deltaPhiDeg << endl;
    // if (entry==10) break;
    
     hPhi_time->Fill(header->realTime, deltaPhiDeg);
     hTheta_time->Fill(header->realTime, deltaThetaDeg);

     
  }


  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1");
  hPhi_time->SetTitle(Form("#phi_{expected}-#phi_{zoom}: %s pulses;Time;#delta #phi (degrees)", arSample.c_str()));
  hPhi_time->GetXaxis()->SetTimeDisplay(1);
  //  hPhi_time->GetXaxis()->SetTimeFormat("%b %d, %H:%M");
  hPhi_time->GetXaxis()->SetTimeFormat("%H:%M");
  //  hPhi_time->GetXaxis()->SetNdivisions(22,12,0, true);
  hPhi_time->Draw("colz");
  
  c1->Print(Form("AngResVStime_phi_%ssample_%sphasecentres.png", arSample.c_str(), pcSample.c_str()));
  c1->Print(Form("AngResVStime_phi_%ssample_%sphasecentres.pdf", arSample.c_str(), pcSample.c_str()));


  hTheta_time->SetTitle(Form("#theta_{expected}-#theta_{zoom}: %s pulses;Time;#delta #theta (degrees)", arSample.c_str()));
  hTheta_time->GetXaxis()->SetTimeDisplay(1);
  //  hTheta_time->GetXaxis()->SetTimeFormat("%b %d, %H:%M");
  hTheta_time->GetXaxis()->SetTimeFormat("%H:%M");
  //  hTheta_time->GetXaxis()->SetNdivisions(22,12,0, true);
  hTheta_time->Draw("colz");
  
  c1->Print(Form("AngResVStime_theta_%ssample_%sphasecentres.png", arSample.c_str(), pcSample.c_str()));
  c1->Print(Form("AngResVStime_theta_%ssample_%sphasecentres.pdf", arSample.c_str(), pcSample.c_str()));


  

}
