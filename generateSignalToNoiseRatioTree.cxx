// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Use reconstructed WAIS arrival direction to coherently sum waveform to find SNR.
********************************************************************************************************* */

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

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

  // Photogrammetry positions
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(Int_t surf=0; surf<NUM_SURF; surf++){
    for(Int_t chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
    }
  }
  TString lindaFileName = "photogrammetryNoExtraDelay";
  

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  TChain* angResChain = new TChain("angResTree");

  OutputConvention oc(argc, argv);
  TString outputDir = oc.getOutputDir();
  if(outputDir.Contains("filterOn")){
    CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite Notch", 260 - 26, 260 + 26);
    CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch", 370 - 26, 370 + 26);
    cc->addNotch(notch260);
    cc->addNotch(notch370);
  }
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
    
    fileName = oc.getOutputDir() + TString::Format("generateAngularResolutionTreePlots_%d*.root", run);
    // fileName = TString::Format("generateAngularResolutionTreePlots_%d*.root", run);        
    angResChain->Add(fileName);
  }
  
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  headChain->BuildIndex("eventNumber");
  gpsChain->BuildIndex("eventNumber");
  calEventChain->BuildIndex("eventNumber");
  
  Double_t peakThetaDeg, peakPhiDeg;
  UInt_t eventNumber;
  // angResChain->SetBranchAddress("zoomPhiDeg3", &peakPhiDeg);
  // angResChain->SetBranchAddress("zoomThetaDeg3", &peakThetaDeg);
  // angResChain->SetBranchAddress("eventNumber", &eventNumber);
  angResChain->SetBranchAddress("zoomPhiDeg", &peakPhiDeg);
  angResChain->SetBranchAddress("zoomThetaDeg", &peakThetaDeg);
  angResChain->SetBranchAddress("eventNumber", &eventNumber);

  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  

  lindaFileNameReference->Write();
  

  TTree* snrTree = new TTree("snrTree", "snrTree");
  Double_t snr2 = 0;
  Double_t snr1 = 0;
  Double_t snr0 = 0;
  snrTree->Branch("snr0", &snr0);
  snrTree->Branch("snr1", &snr1);
  snrTree->Branch("snr2", &snr2);
  UInt_t eventNumberOutput = 0;
  snrTree->Branch("eventNumber", &eventNumberOutput);    
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  

  // lindaFileNameReference->Write();

  Long64_t nEntries = angResChain->GetEntries();
  Long64_t maxEntry = 0; //5; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 5;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    angResChain->GetEntry(entry);
    // std::cout << angResChain->GetEntries() << std::endl;
    // std::cout << eventNumber << std::endl;
    // // eventNumber = 9390257;
    // eventNumber = 60832108;
    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(eventNumber);
    calEventChain->GetEntryWithIndex(eventNumber);

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
    cc->getNormalizedInterpolatedTGraphs(usefulEvent, pol);
    cc->doFFTs(pol);
    Int_t numFreqs = FancyFFTs::getNumFreqs(cc->numSamples);
    Int_t numFreqsPadded = FancyFFTs::getNumFreqs(cc->numSamplesUpsampled);
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      FancyFFTs::zeroPadFFT(cc->ffts[pol][ant], cc->fftsPadded[pol][ant], numFreqs, numFreqsPadded);
    }

    eventNumberOutput = eventNumber;
    TGraph* gr2 = cc->makeUpsampledCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 2, snr2);
    TGraph* gr1 = cc->makeUpsampledCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 1, snr1);
    TGraph* gr0 = cc->makeUpsampledCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 0, snr0);
    // std::cout << snr0 << "\t" << snr1 << "\t" << snr2 << std::endl;

    TGraph* gr2_2 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 2, snr2);
    TGraph* gr1_2 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 1, snr1);
    TGraph* gr0_2 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 0, snr0);
    // std::cout << snr0 << "\t" << snr1 << "\t" << snr2 << std::endl;

    if(numSaved < maxToSave){
      gr2->SetName(TString::Format("%s_2", gr2->GetName()));
      gr1->SetName(TString::Format("%s_1", gr1->GetName()));
      gr0->SetName(TString::Format("%s_0", gr0->GetName()));

      gr2->Write();
      gr1->Write();
      gr0->Write();

      gr2_2->SetName(TString::Format("%s_2", gr2_2->GetName()));
      gr1_2->SetName(TString::Format("%s_1", gr1_2->GetName()));
      gr0_2->SetName(TString::Format("%s_0", gr0_2->GetName()));

      gr2_2->Write();
      gr1_2->Write();
      gr0_2->Write();

      numSaved++;
    }

    delete gr2;
    delete gr1;
    delete gr0;
    delete gr2_2;
    delete gr1_2;
    delete gr0_2;
    

    delete usefulEvent;

    snrTree->Fill();
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();
  return 0;
}
