// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to look at photogrammetry reconstructed positions in januaryTests/photogrammetryPosition
             And apply varying pitch/roll offsets and look at mean and RMS.
*************************************************************************************************************** */

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


// // Copied from UsefulAdu5Pat.h
// // I need to remove the const from this namespace...
// namespace AnitaStaticAdu5Offsets{
//   const Double_t heading = 0;
//   const Double_t pitch = 0;
//   const Double_t roll = 0;
// }


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


  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* angResChain = new TChain("angResTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);    
    gpsChain->Add(fileName);

    // Will wildcard matching work
    // fileName = TString::Format("januaryTests/photogrammetryNumbers/generateAngularResolutionTreePlots_%d*.root", run);
    fileName = TString::Format("januaryTests/photogrammetryNumbersZeroedChannel16BH_cosminV3Trees/generateAngularResolutionTreePlots_%d*.root", run);    
    // fileName = TString::Format("generateAngularResolutionTreePlots_%d_2016-01-29_*.root", run);    
    angResChain->Add(fileName);
  }
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("eventNumber");  
  UInt_t eventNumberGps = 0;
  gpsChain->SetBranchAddress("eventNumber", &eventNumberGps);

  
  UInt_t eventNumber = 0;
  angResChain->SetBranchAddress("eventNumber", &eventNumber);
  Double_t zoomPhiDeg = 0;
  angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
  Double_t zoomThetaDeg = 0;
  angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);








  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  TTree* pitchRollTree = new TTree("pitchRollTree", "pitchRollTree");

  const Int_t numPitches = 11;
  const Double_t startPitch = 0; //-2; //-0.1; //-2;
  const Double_t deltaPitch = 0.4; //0.1; //0.4;
  Double_t pitchOffsets[numPitches] = {0};
  for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
    pitchOffsets[pitchInd] = startPitch + deltaPitch*pitchInd;
  }
  pitchRollTree->Branch(TString::Format("pitchOffsets[%d]", numPitches), pitchOffsets);

  
  const Int_t numRolls = 11;
  const Double_t startRoll = 0; //-2; //-0.5; //-2;
  const Double_t deltaRoll = 0.4; //0.1; //0.4;
  Double_t rollOffsets[numRolls] = {0};
  for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
    rollOffsets[rollInd] = startRoll + deltaRoll*rollInd;
  }  
  pitchRollTree->Branch(TString::Format("rollOffsets[%d]", numRolls), rollOffsets);  




  Double_t deltaPhiDegs[numPitches][numRolls];
  pitchRollTree->Branch("deltaPhiDegs",deltaPhiDegs,
			TString::Format("deltaPhiDegs[%d][%d]/D", numPitches, numRolls));
  Double_t deltaThetaDegs[numPitches][numRolls];
  pitchRollTree->Branch("deltaThetaDegs",deltaThetaDegs,
			TString::Format("deltaThetaDegs[%d][%d]/D", numPitches, numRolls));


  Double_t thetaExpected[numPitches][numRolls];
  pitchRollTree->Branch("thetaExpected",thetaExpected,
			TString::Format("thetaExpected[%d][%d]/D", numPitches, numRolls));
  Double_t phiExpected[numPitches][numRolls];
  pitchRollTree->Branch("phiExpected",phiExpected,
			TString::Format("phiExpected[%d][%d]/D", numPitches, numRolls));
  UInt_t eventNumberPitchRoll = 0;
  pitchRollTree->Branch("eventNumber", &eventNumberPitchRoll);
  Double_t zoomPhiDegPitchRoll = 0;
  pitchRollTree->Branch("zoomPhiDeg", &zoomPhiDegPitchRoll);
  Double_t zoomThetaDegPitchRoll = 0;
  pitchRollTree->Branch("zoomThetaDeg", &zoomThetaDegPitchRoll);
  
  Long64_t nEntries = angResChain->GetEntries();
  Long64_t maxEntry = 0; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    angResChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(eventNumber);

    // std::cout << eventNumber << "\t" << eventNumberGps << std::endl;
    if(eventNumber != eventNumberGps){
      std::cerr << "Warning! Unable to find eventNumber " << eventNumber << " in gpsChain!" << std::endl;
      continue;
    }

    UsefulAdu5Pat usefulPat(pat);

    eventNumberPitchRoll = eventNumber;
    zoomThetaDegPitchRoll = zoomThetaDeg;
    zoomPhiDegPitchRoll = zoomPhiDeg;    
    
    // for(Int_t headingInd=0; headInd < numHeadings; headingInd++){
    //   usefulPat.heading = += headingOffsets[headingInd];
    //   if(usefulPat.heading < 0){
    // 	usefulPat.heading+=360;
    //   }
    //   else if(usefulPat.heading >= 360){
    // 	usefulPat.heading-=360;
    //   }

      for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){	
	usefulPat.pitch = pitchOffsets[pitchInd];
	for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
	  usefulPat.roll = rollOffsets[rollInd];
	
	  usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected[pitchInd][rollInd], phiExpected[pitchInd][rollInd]);

	  // My silly convention that I've stuck to so far...
	  phiExpected[pitchInd][rollInd]*=TMath::RadToDeg();
	  thetaExpected[pitchInd][rollInd]*=-1*TMath::RadToDeg();

	
	  deltaPhiDegs[pitchInd][rollInd] = RootTools::getDeltaAngleDeg(phiExpected[pitchInd][rollInd], zoomPhiDeg);
	  deltaThetaDegs[pitchInd][rollInd] = RootTools::getDeltaAngleDeg(thetaExpected[pitchInd][rollInd], zoomThetaDeg);

	  // std::cout << usefulPat.pitch << "\t" << usefulPat.roll << "\t" << thetaExpected[pitchInd][rollInd] << "\t" << phiExpected[pitchInd][rollInd] << "\t" << deltaPhiDegs[pitchInd][rollInd] << "\t" << deltaThetaDegs[pitchInd][rollInd] << std::endl;


	
	
	}
      }
    // }

    pitchRollTree->Fill();
    
    p++;
  }

  outFile->Write();
  outFile->Close();
  
  return 0;
}
