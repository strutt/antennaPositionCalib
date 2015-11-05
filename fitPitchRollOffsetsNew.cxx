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
#include <TProfile2D.h>
#include <THnSparse.h>

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
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnitaIIINumbers(1);  

  TChain* deltaTChain = new TChain("deltaTTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    deltaTChain->Add(TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run));
  }
  UInt_t eventNumber = 0;
  Double_t correlationDeltaTs[NUM_COMBOS] = {0};
  Double_t correlationValues[NUM_COMBOS] = {0};  
  Double_t correlationDeltaTsClose[NUM_COMBOS] = {0};
  Double_t correlationValuesClose[NUM_COMBOS] = {0};  
  
  deltaTChain->SetBranchAddress("eventNumber", &eventNumber);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTs[%d]", NUM_COMBOS), correlationDeltaTs);
  deltaTChain->SetBranchAddress(TString::Format("correlationValues[%d]", NUM_COMBOS), correlationValuesClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTsClose[%d]", NUM_COMBOS), correlationDeltaTsClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationValuesClose[%d]", NUM_COMBOS), correlationValues);
  
  Adu5Pat* pat =0;
  TChain* gpsChain = RootTools::getAdu5PatChain(firstRun, lastRun, pat);
  gpsChain->BuildIndex("realTime");

  RawAnitaHeader* header = NULL;
  TChain* headChain = RootTools::getHeadChain(firstRun, lastRun, header);
  headChain->BuildIndex("eventNumber");
  
  
  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");

  Long64_t nEntries = deltaTChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  CrossCorrelator* cc = new CrossCorrelator();
  std::vector<Int_t> combos;
  std::vector<Int_t> ant1s;
  std::vector<Int_t> ant2s;  
  
  for(Int_t combo=0; combo < NUM_COMBOS; combo++){
    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);

    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
  }
  delete cc;
  
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry);

    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    for(int combo=0; combo<NUM_COMBOS; combo++){
      Int_t ant1=ant1s.at(combo);
      Int_t ant2=ant2s.at(combo);

      if(correlationDeltaTsClose[combo]==correlationDeltaTs[combo]){
	UsefulAdu5Pat usefulPat(pat);

	Double_t thetaExpected, phiExpected;
	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected *= TMath::RadToDeg();

	Double_t phi1 = geom->getAntPhiPositionRelToAftFore(ant1)*TMath::RadToDeg();
	Double_t phi2 = geom->getAntPhiPositionRelToAftFore(ant2)*TMath::RadToDeg();
	Double_t deltaAngleDeg1 = RootTools::getDeltaAngleDeg(phi1, phiExpected);
	Double_t deltaAngleDeg2 = RootTools::getDeltaAngleDeg(phi2, phiExpected);

	// std::cout << ant1 << "\t" << phiExpected << "\t" << phi1 << "\t" << deltaAngleDeg1 << std::endl;
	// std::cout << ant2 << "\t" << phiExpected << "\t" << phi2 << "\t" << deltaAngleDeg2 << std::endl;	

	if(TMath::Abs(deltaAngleDeg1) < 22.5 && TMath::Abs(deltaAngleDeg2) < 22.5){
	
	  Double_t dt_e = usefulPat.getDeltaTExpected(ant2, ant1,
						      AnitaLocations::LONGITUDE_WAIS,
						      AnitaLocations::LATITUDE_WAIS,
						      AnitaLocations::ALTITUDE_WAIS);

	  std::cout << correlationDeltaTs[combo] << "\t" << dt_e << std::endl;

	}	
      }
      
    }
    p++;
  }

  outFile->Write();
  outFile->Close();
  return 0;

}
