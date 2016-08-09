#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "TStyle.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"  
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

void stupidMacro(){
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corName[FILENAME_MAX];

  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  Adu5Pat *pat2 =0;
  Adu5Pat *pat3 =0;
  // CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *gpsChain2 = new TChain("adu5PatTree");
  TChain *gpsChain3 = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  // TChain *corChain = new TChain("corTree");

  for (unsigned int run=331;run<355;++run){
      
    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/timedHeadFile%d.root",run,run);
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsEvent%d.root",run,run);
    sprintf(corName, "/unix/anita3/linda/corTrees/corRun_NEW10_HPOL_%d.root", run); 

    headChain->Add(headerName);
    gpsChain->Add(gpsName);
    gpsChain2->Add(gpsName);
    gpsChain3->Add(gpsName);
    // corChain->Add(corName);
    
  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);
  gpsChain2->SetBranchAddress("pat",&pat2);
  gpsChain3->SetBranchAddress("pat",&pat3);
  // corChain->SetBranchAddress("cor",&cor);

  UInt_t gps_eventNumber;
  gpsChain->SetBranchAddress("eventNumber", &gps_eventNumber);
  UInt_t gps2_eventNumber;
  gpsChain2->SetBranchAddress("eventNumber", &gps2_eventNumber);
  UInt_t gps3_eventNumber;
  gpsChain3->SetBranchAddress("eventNumber", &gps3_eventNumber);

  //  headChain->BuildIndex("header->eventNumber");
  //  gpsChain->BuildIndex("eventNumber");
  gpsChain2->BuildIndex("pat->realTime");
  gpsChain3->BuildIndex("pat->realTime");

  // int maxEntry=corChain->GetEntries();
  
  // AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  Double_t phiWave, thetaWave, lower, upper;
  Double_t maxCorrTime, deltaTExpected;
  int ant1, ant2;

  int countIndex = 0;
  int countNotsaved = 0;

  bool save = false;

  headChain->GetEntry(0);
  double firstTS = header->triggerTime;
  headChain->GetEntry(headChain->GetEntries()-1);
  double lastTS = header->triggerTime;

  Double_t additionalPhi = 22.5*TMath::DegToRad();
  Double_t TwoPi = TMath::Pi()*2.;
  //  TH2D *h1 = new TH2D("h1", "", 100, firstTS, lastTS, 100, 0, 3600*24);
  TH2D *h1 = new TH2D("h1", "", 1000, firstTS, lastTS, 100, -2, 2);
  TH2D *h2 = new TH2D("h2", "", 1000, firstTS, lastTS, 100, -5, 5);
  TH1D *dHeading1 = new TH1D("dHeading1", "", 100, -2, 2);
  TH1D *dHeading2 = new TH1D("dHeading2", "", 100, -5, 5);
  TH1D *dtrigger = new TH1D("trigger", "", 100, -2e3, 2e3);
  const Double_t maxDeltaTriggerTimeNs = 1200;

  TH1D *dPhi1 = new TH1D("dPhi1", "", 100, -5, 5);
  TH1D *dPhi2 = new TH1D("dPhi2", "", 100, -5, 5);
  TH1D *dTheta1 = new TH1D("dTheta1", "", 100, -2, 2);
  TH1D *dTheta2 = new TH1D("dTheta2", "", 100, -2, 2);



  cout << headChain->GetEntries() << " " << gpsChain->GetEntries() << " " << gpsChain2->GetEntries() << endl;

  double phiWave1, thetaWave1;
  double phiWave2, thetaWave2;
  double phiWave3, thetaWave3;

  for(Long64_t entry=0;entry<headChain->GetEntries();++entry) {
  //  for(Long64_t entry=0;entry<5;++entry) {

    // save=false;
    // corChain->GetEntry(entry);

    // Long64_t headEntry = headChain->GetEntryNumberWithIndex(cor->eventNumber);
    // if(headEntry < 0 ) continue;
   headChain->GetEntry(entry);    
   if((header->trigType & 1)!=1) continue;
   // Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->eventNumber);
   // if(gpsEntry < 0 ) continue;
   gpsChain->GetEntry(entry);

   UsefulAdu5Pat usefulPat(pat);
   int triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
   int triggerTimeNs = header->triggerTimeNs;
   Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
   dtrigger->Fill(deltaTriggerTimeNs);

   if(TMath::Abs(deltaTriggerTimeNs) > maxDeltaTriggerTimeNs) continue;

   usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave1, phiWave1);
    
   Long64_t gpsEntry2 = gpsChain2->GetEntryNumberWithIndex(header->triggerTime);
   if(gpsEntry2 < 0 ) continue;
   gpsChain2->GetEntry(gpsEntry2);

   h1->Fill(header->triggerTime, (pat->heading-pat2->heading) );
   dHeading1->Fill((pat->heading-pat2->heading));
   // if (header->eventNumber!=gps_eventNumber) cout << header->eventNumber << " 1 " << gps_eventNumber << endl; 
   // if (header->eventNumber!=gps2_eventNumber) cout << header->eventNumber << " 2 " << gps2_eventNumber << endl; 

   UsefulAdu5Pat usefulPat2(pat2);
   usefulPat2.getThetaAndPhiWaveWaisDivide(thetaWave2, phiWave2);

   dPhi1->Fill((phiWave1-phiWave2)*TMath::RadToDeg());
   dTheta1->Fill((thetaWave1-thetaWave2)*TMath::RadToDeg());

   Long64_t gpsEntry3 = gpsChain3->GetEntryNumberWithIndex(header->realTime);
   if(gpsEntry3 < 0 ) continue;
   gpsChain3->GetEntry(gpsEntry3);


   UsefulAdu5Pat usefulPat3(pat3);
   usefulPat3.getThetaAndPhiWaveWaisDivide(thetaWave3, phiWave3);


   //    if ((header->triggerTime%(3600*24))!=(pat->timeOfDay/1e3)) cout << header->triggerTime%(3600*24) << " " << pat->timeOfDay/1e3 << endl;
   // double tday = (header->triggerTime%(3600*24));
   // h1->Fill(header->triggerTime, (tday - pat->timeOfDay/1e3) );

   h2->Fill(header->triggerTime, (pat->heading-pat3->heading) );

   // cout << pat->heading << " " << pat2->heading << endl;
    
   dHeading2->Fill((pat->heading-pat3->heading));

   dPhi2->Fill((phiWave1-phiWave3)*TMath::RadToDeg());
   dTheta2->Fill((thetaWave1-thetaWave3)*TMath::RadToDeg());

  }

  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(1);
  dHeading1->SetTitle("Heading(eventNumber)-heading(triggerTime);Heading(eventNumber)-heading(triggerTime) [degrees];Number of events");
  dHeading1->Draw();
  c1->Print("dheading_triggerTime.png");
  dHeading2->SetTitle("Heading(eventNumber)-heading(realTime);Heading(eventNumber)-heading(realTime) [degrees];Number of events");
  dHeading2->Draw();
  c1->Print("dheading_realTime.png");

  h1->SetTitle("Heading(eventNumber)-heading(triggerTime);triggerTime;Heading(eventNumber)-heading(triggerTime) [degrees]");
  h1->Draw("colz");
  c1->Print("dheading_triggerTime_time.png");

  h2->SetTitle("Heading(eventNumber)-heading(realTime);triggerTime;Heading(eventNumber)-heading(realTime) [degrees]");
  h2->Draw("colz");
  c1->Print("dheading_realTime_time.png");


  c1->SetLogy();
  dtrigger->Draw();
  c1->Print("TriggerTimeNsDiff.png");

  dPhi1->SetTitle("PhiExpected(eventNumber)-phiExpected(triggerTime);PhiExpected(eventNumber)-phiExpected(triggerTime) [degrees];Events");
  dPhi1->Draw();
  c1->Print("DeltaPhiExpected_triggerTime.png");


  dPhi2->SetTitle("PhiExpected(eventNumber)-phiExpected(realTime);PhiExpected(eventNumber)-phiExpected(realTime) [degrees];Events");
  dPhi2->Draw();
  c1->Print("DeltaPhiExpected_realTime.png");


  dTheta1->SetTitle("ThetaExpected(eventNumber)-thetaExpected(triggerTime);ThetaExpected(eventNumber)-thetaExpected(triggerTime) [degrees];Events");
  dTheta1->Draw();
  c1->Print("DeltaThetaExpected_triggerTime.png");


  dTheta2->SetTitle("ThetaExpected(eventNumber)-thetaExpected(realTime);ThetaExpected(eventNumber)-thetaExpected(realTime) [degrees];Events");
  dTheta2->Draw();
  c1->Print("DeltaThetaExpected_realTime.png");

}
