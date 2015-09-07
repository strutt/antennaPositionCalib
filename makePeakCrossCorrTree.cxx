#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>
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

  if(argc!=2){
    std::cerr << "Usage:" << argv[0] << "[run]" << std::endl;
  }
  const Int_t run = atoi(argv[1]);

  const Double_t sourceLat = - (79 + (27.93728/60));
  const Double_t sourceLon = -(112 + (6.74974/60));
  const Double_t sourceAlt = 1813.42;

  const Double_t cutTimeNs = 1200; //60e3;

  // TString eventFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/eventFile%d.root", run, run);
  TString eventFileName = TString::Format("~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", 
					  run, run);
  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  TString headFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
  TFile* headFile = TFile::Open(headFileName);
  TTree* headTree = (TTree*) headFile->Get("headTree");

  TString gpsFileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
  TFile* gpsFile = TFile::Open(gpsFileName);
  TTree* adu5PatTree = (TTree*) gpsFile->Get("adu5PatTree");

  // RawAnitaEvent* event = NULL;
  CalibratedAnitaEvent* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  adu5PatTree->SetBranchAddress("pat", &pat);

  FancyTTreeInterpolator gpsInterp(adu5PatTree, "pat->realTime");
  gpsInterp.add("pat->heading", "pat->heading>-500", 360.0);
  gpsInterp.add("pat->altitude");
  gpsInterp.add("pat->latitude");
  gpsInterp.add("pat->longitude");

  TFile* outFile = new TFile(TString::Format("%sPlots%d.root", argv[0], run),"recreate");
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
  Double_t heading = 0;
  Double_t thetaDegExpected = 0;
  Double_t phiDegExpected = 0;

  antTree->Branch("eventNumber", &eventNumber);
  antTree->Branch("peakCorr", &peakCorr);
  antTree->Branch("peakTime", &peakTime);
  antTree->Branch("heading", &heading);
  antTree->Branch("thetaDegExpected", &thetaDegExpected);
  antTree->Branch("phiDegExpected", &phiDegExpected);

  Long64_t nEntries = headTree->GetEntries();
  Long64_t maxEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << std::endl;

  Int_t ant1 = 0;
  Int_t ant2 = 18;

  ProgressBar p(maxEntry);

  TString histName = TString::Format("hTrigTime%d", run);
  TString histTitle = TString::Format("Trigger time (ns), run %d; eventNumber; triggerTime (ns)", run);
  headTree->GetEntry(0);
  TH2D* hTrigTime = new TH2D(histName, histTitle, 1024, 
                             header->eventNumber, maxEntry+header->eventNumber, 
                             10000, 0, 1e9);

  TString histName2 = TString::Format("hTrigTimeZoom%d", run);
  TH2D* hTrigTime2 = new TH2D(histName2, histTitle, 1024, 
                              header->eventNumber, maxEntry+header->eventNumber, 
                              1024, 0, 2e6);



  for(Long64_t entry=0; entry<maxEntry; entry++){

    headTree->GetEntry(entry);
    if((header->trigType & 1)==1){ // RF trigger
      
      Int_t triggerTimeNs = header->triggerTimeNs;
      eventNumber = header->eventNumber;

      hTrigTime->Fill(eventNumber, triggerTimeNs);

      Adu5Pat pat2;
      pat2.heading = gpsInterp.interp("pat->heading", header->realTime);
      pat2.latitude = gpsInterp.interp("pat->latitude", header->realTime);
      pat2.longitude = gpsInterp.interp("pat->longitude", header->realTime);
      pat2.altitude = gpsInterp.interp("pat->altitude", header->realTime);

      UsefulAdu5Pat usefulPat(&pat2);
      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

      Double_t thetaWave, phiWave;
      usefulPat.getThetaAndPhiWaveAnita3(sourceLon, sourceLat, sourceAlt, thetaWave, phiWave);
      thetaDegExpected = thetaWave*TMath::RadToDeg();
      phiDegExpected = phiWave*TMath::RadToDeg();

      if(TMath::Abs(Double_t(triggerTimeNsExpected) - Double_t(triggerTimeNs)) < cutTimeNs){

	eventTree->GetEntry(entry);
	hTrigTime2->Fill(eventNumber, triggerTimeNs);
	heading = pat2.heading;

	// UsefulAnitaEvent usefulEvent(event, WaveCalType::kDefault, header);
	UsefulAnitaEvent usefulEvent(event, WaveCalType::kDefault);//, header);

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
