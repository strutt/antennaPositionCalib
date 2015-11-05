// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to find peak cross correlation offsets between antenna pairs in pulses from Wais Divide.
*************************************************************************************************************** */

#include <TFile.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>


#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

int main(int argc, char *argv[])
{

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;    
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1] << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;
  const Double_t maxDeltaTriggerTimeNs = 1200;


  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
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
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);


  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  Int_t useKurt = 0;
  geom->useKurtAnitaIIINumbers(useKurt); //1;
  

  TString titleSuffix = useKurt==1 ? "Photo" : "Feed";
  TString outFileName = TString::Format("%s_run%d-%d", argv[0], firstRun, lastRun) + titleSuffix + "Plots.root";
  TFile* outFile = new TFile(outFileName, "recreate");    


  const Double_t maxDistKm = 2e3;
  const Int_t numDistBins = 2048;

  const Double_t minThetaDeg = -50;
  const Double_t maxThetaDeg = 50;
  const Int_t numThetaBins = 512; //1024;

  const Double_t minPhiDeg = 0;
  const Double_t maxPhiDeg = 360;
  const Int_t numPhiBins = 1024;
  
  TString title1 = TString::Format("Expected elevation vs distance from WAIS divide runs %d - %d; Distance (km); Elevation (Degrees)", firstRun, lastRun);
  TH2D* hDistEl = new TH2D("hDistEl", title1,
			   numDistBins, 0, maxDistKm, numThetaBins, minThetaDeg, maxThetaDeg);


  TString title2 = TString::Format("WAIS pulses expected arrival directions runs %d - %d; Azimuth (Degrees); Elevation (Degrees)", firstRun, lastRun);
  TH2D* hPhiThetaExpected = new TH2D("hThetaPhiExpected", title2,
				     numPhiBins, minPhiDeg, maxPhiDeg,
				     numThetaBins, minThetaDeg, maxThetaDeg);


  TString title3 = TString::Format("Expected elevation vs heading runs %d - %d; Heading (degrees); Elevation (Degrees)", firstRun, lastRun);
  TH2D* hHeadEl = new TH2D("hHeadEl", title3,
			   numPhiBins, minPhiDeg, maxPhiDeg,
			   numThetaBins, minThetaDeg, maxThetaDeg);  

  Int_t distCutLow = 377;
  Int_t distCutHigh = 383;
  
  TString title4 = TString::Format("Expected elevation vs heading runs %d - %d (%d < distance (km) < %d) ; Heading (degrees); Elevation (Degrees)", firstRun, lastRun, distCutLow, distCutHigh);
  TH2D* hHeadEl2 = new TH2D("hHeadEl2", title4,
			   numPhiBins, minPhiDeg, maxPhiDeg,
			   numThetaBins, minThetaDeg, maxThetaDeg);

 

  TString title5 = TString::Format("Heading vs. Pitch, runs %d - %d", firstRun, lastRun);
  title5 += "; Heading (degrees); Pitch (Degrees)";
  TH2D* hHeadPitch = new TH2D("hHeadPitch", title5,
			      numPhiBins, minPhiDeg, maxPhiDeg,
			      numThetaBins, minThetaDeg, maxThetaDeg);

  TString title6 = TString::Format("Heading vs. Roll, runs %d - %d", firstRun, lastRun);
  title6 += "; Heading (degrees); Roll (Degrees)";
  TH2D* hHeadRoll = new TH2D("hHeadRoll", title6,
			      numPhiBins, minPhiDeg, maxPhiDeg,
			      numThetaBins, minThetaDeg, maxThetaDeg);

  TString title7 = TString::Format("Pitch vs. Roll, runs %d - %d", firstRun, lastRun);
  title7 += "; Pitch (degrees); Roll (Degrees)";
  TH2D* hPitchRoll = new TH2D("hPitchRoll", title7,
			      numThetaBins, minThetaDeg, maxThetaDeg,			      
			      numThetaBins, minThetaDeg, maxThetaDeg);
  



  title5 = TString::Format("Heading vs. Pitch, runs %d - %d (%d < distance (km) < %d)", firstRun, lastRun, distCutLow, distCutHigh);
  title5 += "; Heading (degrees); Pitch (Degrees)";
  TH2D* hHeadPitch2 = new TH2D("hHeadPitch2", title5,
			      numPhiBins, minPhiDeg, maxPhiDeg,
			      numPhiBins, minThetaDeg, maxThetaDeg);

  title6 = TString::Format("Heading vs. Roll, runs %d - %d (%d < distance (km) < %d)", firstRun, lastRun, distCutLow, distCutHigh);
  title6 += "; Heading (degrees); Roll (Degrees)";
  TH2D* hHeadRoll2 = new TH2D("hHeadRoll2", title6,
			      numPhiBins, minPhiDeg, maxPhiDeg,
			      numThetaBins, minThetaDeg, maxThetaDeg);

  title7 = TString::Format("Pitch vs. Roll, runs %d - %d (%d < distance (km) < %d)d", firstRun, lastRun, distCutLow, distCutHigh);
  title7 += "; Pitch (degrees); Roll (Degrees)";
  TH2D* hPitchRoll2 = new TH2D("hPitchRoll2", title7,
			      numThetaBins, minThetaDeg, maxThetaDeg,			      
			      numThetaBins, minThetaDeg, maxThetaDeg);

  
  
  UInt_t eventNumber;
  // Double_t heading;
  Double_t thetaExpected;
  Double_t phiExpected;
  UInt_t triggerTimeNs;
  UInt_t triggerTimeNsExpected;    

    
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
	eventNumber = header->eventNumber;
	Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;

	// std::cout << pat->heading << "\t" << pat->pitch << "\t" << pat->roll << std::endl;

	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();


	hHeadPitch->Fill(pat->heading, pat->pitch);
	hHeadRoll->Fill(pat->heading, pat->roll);
	hPitchRoll->Fill(pat->pitch, pat->roll);
	

	
	hDistEl->Fill(distKm, thetaExpected);
	hHeadEl->Fill(usefulPat.heading, thetaExpected);
	hPhiThetaExpected->Fill(phiExpected, thetaExpected);

	if(distKm > distCutLow && distKm < distCutHigh){
	  hHeadEl2->Fill(usefulPat.heading, thetaExpected);
	  hHeadPitch2->Fill(pat->heading, pat->pitch);
	  hHeadRoll2->Fill(pat->heading, pat->roll);
	  hPitchRoll2->Fill(pat->pitch, pat->roll);
	  
	}
	
		
      }
    }
    p++;
  }
  

  outFile->Write();
  outFile->Close();

  return 0;
}
