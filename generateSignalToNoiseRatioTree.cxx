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

  TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";
  Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);  
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  TChain* angResChain = new TChain("angResTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/antennaPositionCalib/januaryTests/testing/comboStudies/generateAngularResolutionTreePlots_%d*.root", run);
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
  angResChain->SetBranchAddress("zoomPhiDeg3", &peakPhiDeg);
  angResChain->SetBranchAddress("zoomThetaDeg3", &peakThetaDeg);
  angResChain->SetBranchAddress("eventNumber", &eventNumber);
    
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  

  lindaFileNameReference->Write();

  angResChain->GetEntry(0);
  // eventNumber = 9390257;
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

  // peakPhiDeg = 245.25;
  // peakThetaDeg = -6.05;
  
  TGraph* gr2 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 2);
  TGraph* gr1 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 1);
  TGraph* gr0 = cc->makeCoherentlySummedWaveform(pol, peakPhiDeg, peakThetaDeg, 0);

  gr2->SetName("gr2");
  gr1->SetName("gr1");
  gr0->SetName("gr0");
  
  gr2->Write();
  gr1->Write();
  gr0->Write();
}
