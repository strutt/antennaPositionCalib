#include <UsefulAnitaEvent.h>
#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <FFTtools.h>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include <iostream>
#include <vector>
#include <assert.h>

#include <RootTools.h>
#include <ProgressBar.h>
#include <FancyTTreeInterpolator.h>

int main(int argc, char* argv[]){

  if(argc>1){
    
  }

  const Int_t startRun = 331;
  const Int_t endRun = 355;

  // Thanks steph
  const Double_t sourceLat = - (79 + (27.93728/60));
  const Double_t sourceLon = -(112 + (6.74974/60));
  const Double_t sourceAlt = 1813.42;

  const Double_t cutTimeNs = 1200; //60e3;

  TChain* headChain = new TChain("headTree");
  TChain* eventChain = new TChain("eventTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  for(Int_t run=startRun; run<endRun; run++){
    TString eventFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/eventFile%d.root", run, run);
    eventChain->Add(eventFileName);

    TString rawHeaderFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/eventHeadFile%d.root", 
						run, run);
    headChain->Add(rawHeaderFileName);

    TString gpsFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(gpsFileName);
  }

  RawAnitaEvent* event = NULL;
  eventChain->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  FancyTTreeInterpolator gpsInterp(gpsChain, "pat->realTime");
  gpsInterp.add("pat->heading", "pat->heading>-500", 360.0);
  gpsInterp.add("pat->altitude");
  gpsInterp.add("pat->latitude");
  gpsInterp.add("pat->longitude");

  TFile* outFile = new TFile(TString::Format("%sPlots.root", argv[0]),"recreate");
  gpsInterp.get("pat->heading")->SetName("grHead");
  gpsInterp.get("pat->altitude")->SetName("grAlt");
  gpsInterp.get("pat->latitude")->SetName("grLat");
  gpsInterp.get("pat->longitude")->SetName("grLong");
  
  // TGraph* gr = gpsInterp.get("pat->heading");
  // gpsInterp.get("pat->heading")->Write();
  // gpsInterp.get("pat->altitude")->Write();
  // gpsInterp.get("pat->latitude")->Write();
  // gpsInterp.get("pat->longitude")->Write();

  TTree* antTree = new TTree("antTree","antTree");
  UInt_t eventNumber = 0;
  Double_t peakCorr = 0;
  Double_t peakTime = 0;
  antTree->Branch("eventNumber", &eventNumber);
  antTree->Branch("peakCorr", &peakCorr);
  antTree->Branch("peakTime", &peakTime);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 10000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << std::endl;

  Int_t ant1 = 0;
  Int_t ant2 = 18;

  ProgressBar p(maxEntry);

  // TH2D* hTrigTime = ;

  for(Long64_t entry=0; entry<maxEntry; entry++){

    headChain->GetEntry(entry);
    
    if((header->trigType & 1)==1){ // RF trigger
      
      Int_t triggerTimeNs = header->triggerTimeNs;
      Adu5Pat pat2;
      pat2.heading = gpsInterp.interp("pat->heading", header->realTime);
      pat2.latitude = gpsInterp.interp("pat->latitude", header->realTime);
      pat2.longitude = gpsInterp.interp("pat->longitude", header->realTime);
      pat2.altitude = gpsInterp.interp("pat->altitude", header->realTime);

      UsefulAdu5Pat usefulPat(&pat2);
      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

      if(TMath::Abs(Double_t(triggerTimeNsExpected) - Double_t(triggerTimeNs)) < cutTimeNs){

	eventChain->GetEntry(entry);
      
	UsefulAnitaEvent usefulEvent(event, WaveCalType::kDefault, header);
      

	TGraph* gr1 = usefulEvent.getGraph(ant1);
	TGraph* gr2 = usefulEvent.getGraph(ant2);

	TGraph* grCorr = FFTtools::getInterpolatedCorrelationGraph(gr1, gr2, 0.01);

	Double_t minX, minY;
	RootTools::getMaxMin(grCorr, peakCorr, peakTime,minY, minX);

	antTree->Fill();

	delete gr1;
	delete gr2;
	delete grCorr;
      }
    }  
    p++;
  }
  
  outFile->Write();
  outFile->Close();
  
}
