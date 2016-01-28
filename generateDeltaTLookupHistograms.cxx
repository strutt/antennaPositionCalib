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

  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1] << "\t" << argv[2] << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);

  TChain* deltaTChain = new TChain("deltaTTree");

  Adu5Pat* pat = NULL;
  TChain* gpsChain = RootTools::getAdu5PatChain(firstRun, lastRun, pat);

  RawAnitaHeader* header = NULL;
  TChain* headChain = RootTools::getHeadChain(firstRun, lastRun, header);

  headChain->BuildIndex("eventNumber");
  gpsChain->BuildIndex("realTime");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run);
    deltaTChain->Add(fileName);
  }

  UInt_t eventNumber = 0;
  Double_t correlationDeltaTs[NUM_COMBOS] = {0};
  Double_t correlationValues[NUM_COMBOS] = {0};  
  Double_t correlationDeltaTsClose[NUM_COMBOS] = {0};
  Double_t correlationValuesClose[NUM_COMBOS] = {0};  
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  deltaTChain->SetBranchAddress("eventNumber", &eventNumber);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTs[%d]", NUM_COMBOS), correlationDeltaTs);
  deltaTChain->SetBranchAddress(TString::Format("correlationValues[%d]", NUM_COMBOS), correlationValuesClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTsClose[%d]", NUM_COMBOS), correlationDeltaTsClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationValuesClose[%d]", NUM_COMBOS), correlationValues);
  deltaTChain->SetBranchAddress("thetaExpected", &thetaExpected);
  deltaTChain->SetBranchAddress("phiExpected", &phiExpected);
  
  Long64_t nEntries = deltaTChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");

  const Double_t phiDegMin = 0;
  const Double_t phiDegMax = 360;
  const Int_t numBinsPhi = 2048; //3600
  const Double_t thetaDegMin = -20;
  const Double_t thetaDegMax = 0;
  const Int_t numBinsTheta = 128; //200;
  const Double_t maxDeltaPhiDeg = 22.5; //2;

  const Double_t minDeltaTBin = -50;
  const Double_t maxDeltaTBin = 50;
  const Int_t numDeltaTBins = 1024;

  const Double_t minCorrBin = 0;
  const Double_t maxCorrBin = 1;
  const Int_t numCorrBins = 1024;

  const Int_t numCombos = NUM_COMBOS;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  
  CrossCorrelator* cc = new CrossCorrelator();
  std::vector<Int_t> combos;
  std::vector<Int_t> ant1s;
  std::vector<Int_t> ant2s;  
  TProfile2D* hDtProfs[numCombos];
  // TH2D* hThetaPhiExpecteds[numCombos];
  // TH2D* hPhiExpPhiAnts[numCombos];
  // TH2D* hCorrDts[numCombos];

  THnSparseF* hThetaPhiExpecteds[numCombos];
  THnSparseF* hPhiExpPhiAnts[numCombos];
  THnSparseF* hCorrDts[numCombos];

  TGraph* gr[NUM_COMBOS];

  for(Int_t combo=0; combo < numCombos; combo++){

    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);

    gr[combo] = new TGraph();
    gr[combo]->SetName(TString::Format("grCombo%d", combo));
    

    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
    
    TString name = TString::Format("hDtProf_%d_%d", ant1, ant2);
    TString title = TString::Format("2D profile of #deltat_{measured} as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);
    hDtProfs[combo] = new TProfile2D(name, title,
				     numBinsPhi, phiDegMin, phiDegMax,
				     numBinsTheta, thetaDegMin, thetaDegMax);

    name = TString::Format("hThetaPhiExpected_%d_%d", ant1, ant2);
    title = TString::Format("Number of events as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);

    const Int_t nDim1 = 2;
    Int_t nBins1[nDim1] = {numBinsPhi, numBinsTheta};
    Double_t mins1[nDim1] = {phiDegMin, thetaDegMin};
    Double_t maxs1[nDim1] = {phiDegMax, thetaDegMax};
    
    hThetaPhiExpecteds[combo] = new THnSparseF(name,title, nDim1, nBins1, mins1, maxs1);

    name = TString::Format("hPhiExpPhiAnt_%d_%d", ant1, ant2);
    title = TString::Format("Antenna azimuth vs. expected azimuth for antennas %d and %d; Expected azimuth (degrees); Antenna azimuth (degrees)", ant1, ant2);
    const Int_t nDim2 = 2;
    Int_t nBins2[nDim2] = {numBinsPhi, numBinsTheta};
    Double_t mins2[nDim2] = {phiDegMin, phiDegMin};
    Double_t maxs2[nDim2] = {phiDegMax, phiDegMax};

    hPhiExpPhiAnts[combo] = new THnSparseF(name, title, nDim2, nBins2, mins2, maxs2);

    name = TString::Format("hCorrDts_%d_%d", ant1, ant2);
    title = TString::Format("Correlation coefficient vs #deltat_{measured} between antennas %d and %d; #deltat_{measured} (ns); Correlation coefficient (no units)", ant1, ant2);
    const Int_t nDim3 = 2;
    Int_t nBins3[nDim3] = {numDeltaTBins, numCorrBins};
    Double_t mins3[nDim3] = {minDeltaTBin, minCorrBin};
    Double_t maxs3[nDim3] = {maxDeltaTBin, maxCorrBin};
    hCorrDts[combo] = new THnSparseF(name, title, nDim3, nBins3, mins3, maxs3);
    
  }
  delete cc;
  
  std::cout << "Generating lookup TProfile2Ds" << std::endl;
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry);

    headChain->GetEntryWithIndex(eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);

    Double_t phiExpected0, thetaExpected0;
    UsefulAdu5Pat usefulPat(pat);
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0, phiExpected0);
    // usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0, phiExpected0);
    thetaExpected0 *= -1*TMath::RadToDeg();
    phiExpected0 *= TMath::RadToDeg();

    for(Int_t combo=0; combo < numCombos; combo++){
      Int_t ant1 = ant1s.at(combo);
      Int_t ant2 = ant2s.at(combo);

      Double_t antPhiDeg1 = geom->getAntPhiPositionRelToAftFore(ant1)*TMath::RadToDeg();
      Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(phiExpected0,
							  antPhiDeg1);

      Double_t antPhiDeg2 = geom->getAntPhiPositionRelToAftFore(ant2)*TMath::RadToDeg();
      Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiExpected0,
							  antPhiDeg2);
      
      if(TMath::Abs(deltaPhiDeg1) < maxDeltaPhiDeg && TMath::Abs(deltaPhiDeg2) < maxDeltaPhiDeg){
	
	if(correlationDeltaTsClose[combo]==correlationDeltaTs[combo]){

	  Double_t dt_m = correlationDeltaTsClose[combo];
	  Double_t dt_e = usefulPat.getDeltaTExpected(ant2, ant1,
	  					      phiExpected0*TMath::DegToRad(),
	  					      -1*thetaExpected0*TMath::DegToRad());

	  gr[combo]->SetPoint(gr[combo]->GetN(), phiExpected, dt_e - dt_m);

	  hDtProfs[combo]->Fill(phiExpected0, thetaExpected0, correlationDeltaTsClose[combo]);

	  const Int_t nDim = 2;
	  Double_t coords1[nDim] = {phiExpected0, thetaExpected0};
	  hThetaPhiExpecteds[combo]->Fill(coords1);

	  Double_t coords2[nDim] = {phiExpected0, antPhiDeg1};
	  hPhiExpPhiAnts[combo]->Fill(coords2);

	  Double_t coords3[nDim] = {phiExpected0, antPhiDeg2};
	  hPhiExpPhiAnts[combo]->Fill(coords3);

	  const Int_t nDim4 = 2;
	  Double_t coords4[nDim4] = {correlationDeltaTs[combo], correlationValues[combo]};
	  hCorrDts[combo]->Fill(coords4);
	}
      }
    }
    p++;
  }
  
  THnSparseF* hSparses[numCombos];
  THnSparseF* hSparses2[numCombos];  

  for(Int_t combo=0; combo<numCombos; combo++){
    Int_t ant1 = ant1s.at(combo);
    Int_t ant2 = ant2s.at(combo);

    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    TString title = TString::Format("2D profile of #deltat_{measured} as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);

    const Int_t nDim = 2;
    Int_t nBins[nDim] = {numBinsPhi, numBinsTheta};
    Double_t mins[nDim] = {phiDegMin, thetaDegMin};
    Double_t maxs[nDim] = {phiDegMax, thetaDegMax};
    hSparses[combo] = new THnSparseF(name, title, nDim, nBins, mins, maxs);

    TString name2 = TString::Format("hDtErrorSparse_%d_%d", ant1, ant2);
    TString title2 = TString::Format("2D profile of #sigma#(deltat_{measured}) as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);
    hSparses2[combo] = new THnSparseF(name2, title2, nDim, nBins, mins, maxs);
    // std::cout << name2.Data() << std::endl;

    
    // Double_t halfThetaBinWidth = 0.5*(hDtProfs[combo]->GetYaxis()->GetBinLowEdge(2) - hDtProfs[combo]->GetYaxis()->GetBinLowEdge(1));
    // Double_t halfPhiBinWidth = 0.5*(hDtProfs[combo]->GetXaxis()->GetBinLowEdge(2) - hDtProfs[combo]->GetXaxis()->GetBinLowEdge(1));
    Double_t halfThetaBinWidth = 0;
    Double_t halfPhiBinWidth = 0;

    // std::cout << halfThetaBinWidth << "\t" << halfPhiBinWidth << std::endl;

    for(Int_t binx=1; binx<=hDtProfs[combo]->GetXaxis()->GetNbins(); binx++){
      Double_t phi = hDtProfs[combo]->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
      // phi *= TMath::DegToRad();
      for(Int_t biny=1; biny<=hDtProfs[combo]->GetYaxis()->GetNbins(); biny++){
	if(hDtProfs[combo]->GetBinContent(binx, biny) != 0){
	  Double_t theta = hDtProfs[combo]->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	  // theta *= TMath::DegToRad();
	  Double_t dt_m = hDtProfs[combo]->GetBinContent(binx, biny);
	  Double_t dt_m_error = hDtProfs[combo]->GetBinError(binx, biny);

	  if(dt_m_error > 0){
	    Double_t coords[nDim] = {phi, theta};
	    hSparses[combo]->Fill(coords, dt_m);
	    hSparses2[combo]->Fill(coords, dt_m_error);
	  }
	}
      }
    }
    delete hDtProfs[combo];
    hSparses[combo]->Write();
    hSparses2[combo]->Write();    

    hThetaPhiExpecteds[combo]->Write();
    hPhiExpPhiAnts[combo]->Write();
    hCorrDts[combo]->Write();

    gr[combo]->Write();
    
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
