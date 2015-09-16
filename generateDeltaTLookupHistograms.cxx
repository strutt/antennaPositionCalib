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
  std::cout << argv[0] << "\t" << argv[1] << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);

  TChain* deltaTChain = new TChain("deltaTTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run);
    deltaTChain->Add(fileName);
  }

  std::vector<std::vector<Double_t> >* correlationDeltaTs = NULL;
  std::vector<std::vector<Double_t> >* correlationValues = NULL;  
  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  // std::vector<Double_t>* deltaPhiDeg = NULL;
  deltaTChain->SetBranchAddress("correlationDeltaTs", &correlationDeltaTs);
  deltaTChain->SetBranchAddress("correlationValues", &correlationValues);  
  deltaTChain->SetBranchAddress("thetaExpected", &thetaExpected);
  deltaTChain->SetBranchAddress("phiExpected", &phiExpected);
  // deltaTChain->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);  
  
  Long64_t nEntries = deltaTChain->GetEntries();
  Long64_t maxEntry = 0; //3000;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry);

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  
  // Int_t upsampleFactor = 32;
  // CrossCorrelator* cc = new CrossCorrelator(upsampleFactor);
  // Double_t phiDegMin = 0;
  // Double_t phiDegMax = 360;
  const Double_t phiDegMin = 0;
  const Double_t phiDegMax = 360;
  const Int_t numBinsPhi = 2048;
  const Double_t thetaDegMin = -50;
  const Double_t thetaDegMax = 50;
  const Int_t numBinsTheta = 1024;
  const Double_t maxDeltaPhiDeg = 22.5*2;

  const Double_t minDeltaTBin = -50;
  const Double_t maxDeltaTBin = 50;
  const Int_t numDeltaTBins = 1024;

  const Double_t minCorrBin = 0;
  const Double_t maxCorrBin = 1;
  const Int_t numCorrBins = 1024;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  const Int_t numCombos = NUM_PHI;
  
  CrossCorrelator* cc = new CrossCorrelator();
  std::vector<Int_t> combos;
  std::vector<Int_t> ant1s;
  std::vector<Int_t> ant2s;  
  TProfile2D* hDtProfs[numCombos];
  TH2D* hThetaPhiExpecteds[numCombos];
  TH2D* hPhiExpPhiAnts[numCombos];
  TH2D* hCorrDts[numCombos];

  Double_t lowDtCuts[numCombos] = {2.92969, 3.22266, 2.92969, 3.125,
				   3.22266, 3.32031, 3.02734, 3.22266,
				   3.41797, 3.51562, 2.92969, 3.125,
				   3.02734, 3.22266, 3.41797, 3.22266};

  
  Double_t highDtCuts[numCombos] = {5.46875, 5.85938, 5.37109, 5.85938,
				    5.66406, 5.66406, 5.56641, 5.85938,
				    5.66406, 5.76172, 5.56641, 5.76172,
				    5.37109, 5.95703, 5.76172, 5.76172};  

  
  for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
    Int_t ant1 = comboInd; // i.e. phi
    // Int_t ant2 = (comboInd + 1)%NUM_PHI;
    // Int_t ant2 = comboInd + 2*NUM_PHI;
    Int_t ant2 = comboInd + NUM_PHI;
    Int_t combo = cc->comboIndices[ant1][ant2];
    std:: cout << ant1 << "\t" << ant2 << "\t" << combo << std::endl;
    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
    
    TString name = TString::Format("hDtProf_%d_%d", ant1, ant2);
    TString title = TString::Format("Two dimensional profile of #deltat_{measured} as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);
    hDtProfs[comboInd] = new TProfile2D(name, title,
					numBinsPhi, phiDegMin, phiDegMax,
					numBinsTheta, thetaDegMin, thetaDegMax);

    name = TString::Format("hThetaPhiExpected_%d_%d", ant1, ant2);
    title = TString::Format("Two dimensional histogram counting number of events as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);
    hThetaPhiExpecteds[comboInd] = new TH2D(name,title, 
					    numBinsPhi, phiDegMin, phiDegMax,
					    numBinsTheta, thetaDegMin, thetaDegMax);
    name = TString::Format("hPhiExpPhiAnt_%d_%d", ant1, ant2);
    title = TString::Format("Antenna azimuth vs. expected azimuth for antennas %d and %d; Expected azimuth (degrees); Antenna azimuth (degrees)", ant1, ant2);
    hPhiExpPhiAnts[comboInd] = new TH2D(name, title,
					numBinsPhi, phiDegMin, phiDegMax,
					numBinsPhi, phiDegMin, phiDegMax);
    name = TString::Format("hCorrDts_%d_%d", ant1, ant2); 
    title = TString::Format("Correlation coefficient vs #deltat_{measured} between antennas %d and %d; #deltat_{measured} (ns); Correlation coefficient (no units)", ant1, ant2), 
    hCorrDts[comboInd] = new TH2D(name, title,
				  numDeltaTBins, minDeltaTBin, maxDeltaTBin,
				  numCorrBins, minCorrBin, maxCorrBin);
    
  }
  delete cc;

  
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  std::cout << "Generating lookup TProfile2Ds" << std::endl;
  for(Long64_t entry = 0; entry < maxEntry; entry++){
    deltaTChain->GetEntry(entry);

    for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
      Int_t combo = combos.at(comboInd);
      Int_t ant1 = ant1s.at(comboInd);
      Int_t ant2 = ant2s.at(comboInd);
      
      Double_t antPhiDeg1 = geom->getAntPhiPositionRelToAftFore(ant1)*TMath::RadToDeg();
      Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(phiExpected, antPhiDeg1);

      Double_t antPhiDeg2 = geom->getAntPhiPositionRelToAftFore(ant2)*TMath::RadToDeg();
      Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiExpected, antPhiDeg2);
    
      if(TMath::Abs(deltaPhiDeg1) < maxDeltaPhiDeg && TMath::Abs(deltaPhiDeg2) < maxDeltaPhiDeg){
	if(correlationDeltaTs->at(pol).at(combo) > lowDtCuts[comboInd] && correlationDeltaTs->at(pol).at(combo) < highDtCuts[comboInd]){
	  hDtProfs[comboInd]->Fill(phiExpected, thetaExpected, correlationDeltaTs->at(pol).at(combo));
	  // hDtProf->Fill(deltaPhiDeg->at(0), thetaExpected, correlationDeltaTs->at(pol).at(combo));
	  hThetaPhiExpecteds[comboInd]->Fill(phiExpected, thetaExpected);
	  hPhiExpPhiAnts[comboInd]->Fill(phiExpected, antPhiDeg1);
	  hPhiExpPhiAnts[comboInd]->Fill(phiExpected, antPhiDeg2);	  
	}
	// }
	hCorrDts[comboInd]->Fill(correlationDeltaTs->at(pol).at(combo), correlationValues->at(pol).at(combo));
      }
    }
    p++;
  }

  for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
    TH1D* hCorrDtProj = hCorrDts[comboInd]->ProjectionX();
    Int_t peakBin = RootTools::getPeakBinOfHistogram(hCorrDtProj);
    Double_t peakVal = hCorrDtProj->GetBinContent(peakBin);
    Int_t numBins = hCorrDtProj->GetNbinsX();
    
    // Here I will try finding the first local minima below 1./4 the height of the peak.
    Double_t factor = 0.05;

    Double_t lastVal = peakVal;
    Int_t minBinLow = -1;
    for(int binx=peakBin; binx >=0; binx--){
      Double_t val = hCorrDtProj->GetBinContent(binx);
      if(val < factor*peakVal){
	if(val > lastVal){
	  minBinLow = binx;
	  break;
	}
      }
      lastVal = val;
    }

    lastVal = peakVal;
    Int_t minBinHigh = -1;    
    for(int binx=peakBin; binx <=numBins; binx++){
      Double_t val = hCorrDtProj->GetBinContent(binx);
      if(val < factor*peakVal){
	if(val > lastVal){
	  minBinHigh = binx;
	  break;
	}
      }
      lastVal = val;
    }

    Double_t binWidth = hCorrDtProj->GetBinLowEdge(2) - hCorrDtProj->GetBinLowEdge(1);
    Double_t minVal = hCorrDtProj->GetBinLowEdge(minBinLow) + binWidth;
    Double_t maxVal = hCorrDtProj->GetBinLowEdge(minBinHigh);// - binWidth;

    std::cout << combos.at(comboInd) << "\t"
	      << ant1s.at(comboInd) << "\t"
	      << ant2s.at(comboInd) << "\t"
	      << minVal << "\t" << maxVal
	      << std::endl;
    
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
