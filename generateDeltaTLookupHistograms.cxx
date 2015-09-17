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

  const Int_t numCombos = NUM_COMBOS;
  
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


  Double_t lowDtCuts[numCombos] = {-1.75781, 1.75781, 1.26953, -1.5625, 2.53906, 3.51562, 2.92969, 3.51562, -1.5625, 2.53906, 2.63672, -1.75781, 2.24609, 2.44141, -1.46484, 2.05078, 2.92969, 3.02734, 2.83203, 3.22266, 3.41797, -1.07422, 2.83203, 3.02734, -1.46484, 2.44141, 3.02734, 2.14844, 1.95312, 3.71094, 4.00391, 2.92969, 3.51562, -2.05078, 2.34375, 2.44141, -50.0977, -50.0977, -50.0977, 4.19922, 4.19922, 4.58984, 5.07812, 3.125, 3.61328, -1.26953, 2.63672, 3.125, -1.5625, 2.44141, 2.83203, -50.0977, -50.0977, 2.44141, 2.92969, 3.22266, 3.41797, -1.46484, 2.53906, 2.83203, -1.75781, 2.14844, 2.24609, 1.95312, 0.292969, 2.44141, 2.73438, 3.32031, 3.90625, -1.26953, 2.44141, 0.683594, -4.00391, 0.488281, 2.14844, 2.44141, 2.24609, 2.14844, 2.44141, 3.02734, 3.71094, -1.5625, 2.44141, 2.83203, -1.5625, 2.24609, 2.63672, 2.73438, 2.83203, 2.92969, 3.41797, 3.22266, 4.00391, -1.46484, 2.92969, 3.125, -1.5625, 2.44141, 3.02734, 2.34375, 1.36719, 2.92969, 2.83203, 3.41797, 3.61328, -1.5625, 2.05078, 2.63672, -1.46484, 1.95312, 2.53906, 2.44141, 3.02734, 2.73438, 2.24609, 3.51562, 3.51562, -1.26953, 2.83203, 3.41797, -1.46484, 2.14844, 2.83203, 2.14844, 2.53906, 2.44141, 2.63672, 2.92969, 3.61328, -1.66016, 2.83203, 2.73438, -1.66016, 2.14844, 0.195312, 1.26953, 1.66016, 2.63672, 3.22266, 3.125, 3.61328, -1.07422, 2.73438, 0.585938, -1.26953, 2.14844, 2.53906, 2.24609, 2.63672, 2.73438, 3.125, 3.02734, 3.71094, -1.46484, 2.44141, 2.53906, -1.66016, 2.14844, 2.44141, 2.44141, 2.05078, 2.73438, 0.683594, 3.22266, 3.90625, -1.46484, 2.73438, 2.92969, -1.46484, 2.24609, 3.22266, -0.292969, -2.83203, 1.95312, 1.75781, 3.41797, 4.00391, -1.46484, 2.24609, 2.83203, 1.95312, 0.390625, 2.24609, 2.34375, 3.125, 2.44141, 3.22266, 4.00391, 2.44141, 2.92969, 2.44141, 2.73438, -2.44141, -2.24609, -1.95312, -1.26953, -0.390625, -2.24609, -1.85547, -2.63672, -2.14844, -2.53906, -2.05078, -1.95312, -0.585938, -2.34375, -1.66016, -2.92969, -3.02734, -2.14844, -0.976562, -1.07422, -2.44141, -2.14844, -50.0977, -50.0977, -0.195312, 0, -0.585938, -2.05078, -1.95312, -2.73438, -2.24609, -50.0977, -1.75781, -0.292969, -2.24609, -1.75781, -2.83203, -2.24609, -2.24609, -1.95312, -0.292969, -2.24609, -2.05078, -2.63672, -2.44141, -2.24609, -1.85547, -0.292969, -2.14844, -1.85547, -2.53906, -3.51562, -2.34375, -1.75781, -0.195312, -2.24609, -1.66016, -3.51562, -2.24609, -3.125, -1.85547, -0.585938, -2.44141, -2.05078, -5.17578, -2.24609, -2.24609, -1.85547, -0.488281, -2.05078, -2.05078, -2.63672, -2.14844, -2.44141, -1.66016, -0.878906, -2.34375, -1.75781, -2.73438, -5.27344, -2.24609, -1.95312, -0.488281, -2.34375, -4.78516, -2.83203, -2.44141, -2.24609, -1.85547, -0.292969, -2.63672, -1.46484, -2.63672, -2.44141, -2.34375, -4.49219, -0.878906, -2.24609, -1.85547, -2.73438, -1.95312, -4.88281, -1.75781, -0.488281, -2.34375, -1.36719, -2.63672, -2.34375, -1.95312, -0.292969, -3.22266, -2.44141, -2.73438, -1.75781, -2.34375, -2.83203, -2.83203, -2.05078, -5.17578, -2.44141, -50.0977, -2.14844, -2.53906, -2.24609, -2.63672, -2.34375, -3.02734, -2.14844, -2.63672, -2.14844, -3.41797, -2.24609, -2.73438, -2.05078, -2.92969, -2.63672, -5.56641, -4.78516, -2.63672, -2.14844, -2.63672, -2.24609, -2.34375, -1.85547};
  Double_t highDtCuts[numCombos] = {1.36719, 6.73828, 7.03125, 1.17188, 6.05469, 6.83594, 5.46875, 6.05469, 1.17188, 6.05469, 6.34766, 1.5625, 6.83594, 7.03125, 1.36719, 6.93359, 7.91016, 6.25, 6.73828, 5.85938, 6.25, 1.46484, 5.95703, 6.34766, 0.878906, 5.37109, 5.76172, 6.44531, 7.32422, 6.15234, 6.44531, 5.37109, 6.05469, 0.0976562, 4.88281, 5.17578, -50.1953, -50.1953, -50.1953, 6.54297, 6.93359, 6.15234, 6.64062, 5.85938, 6.25, 1.5625, 6.34766, 6.64062, 3.02734, 7.12891, 7.42188, -50.1953, -50.1953, 6.25, 6.83594, 5.66406, 6.05469, 1.17188, 6.15234, 7.12891, 1.85547, 6.64062, 7.32422, 6.64062, 7.91016, 6.44531, 6.64062, 5.66406, 6.25, 1.5625, 6.54297, 6.73828, 1.5625, 6.73828, 7.32422, 6.83594, 7.42188, 6.25, 6.64062, 5.56641, 6.15234, 1.07422, 6.15234, 6.73828, 1.5625, 6.73828, 7.03125, 7.51953, 8.10547, 6.25, 6.73828, 5.85938, 6.15234, 1.46484, 6.25, 6.83594, 1.36719, 7.12891, 7.03125, 6.83594, 7.22656, 6.34766, 6.93359, 5.66406, 6.05469, 1.26953, 6.15234, 6.54297, 1.75781, 6.73828, 7.03125, 6.83594, 7.22656, 6.34766, 6.93359, 5.76172, 6.15234, 1.5625, 6.34766, 6.64062, 1.46484, 7.03125, 7.22656, 7.61719, 7.12891, 6.25, 6.83594, 5.56641, 6.05469, 1.26953, 6.15234, 6.73828, 1.66016, 6.54297, 7.32422, 7.03125, 7.12891, 6.15234, 6.73828, 5.76172, 6.34766, 1.66016, 6.25, 6.64062, 1.66016, 7.03125, 7.22656, 6.73828, 7.51953, 6.25, 6.44531, 5.37109, 5.85938, 1.36719, 6.44531, 6.54297, 3.90625, 6.83594, 7.42188, 6.54297, 7.03125, 5.95703, 6.64062, 5.95703, 6.05469, 1.46484, 6.25, 6.64062, 1.36719, 6.83594, 7.42188, 7.03125, 9.27734, 6.44531, 6.73828, 5.76172, 6.34766, 1.36719, 6.34766, 7.22656, 6.83594, 7.03125, 6.83594, 7.22656, 6.25, 6.73828, 5.76172, 6.73828, 6.25, 6.73828, 6.93359, 7.03125, 2.73438, 3.22266, 2.24609, 2.83203, 1.36719, 2.34375, 2.63672, 2.83203, 3.02734, 2.92969, 3.61328, 2.83203, 1.5625, 0.976562, 1.46484, 0.488281, 0.878906, 3.71094, 2.53906, 1.36719, 0.292969, 0.390625, -50.1953, -50.1953, 3.125, 3.51562, 1.26953, 2.44141, 2.63672, 2.73438, 3.22266, -50.1953, 2.53906, 1.26953, 1.95312, 2.53906, 2.73438, 2.92969, 3.32031, 2.63672, 1.5625, 2.14844, 2.83203, 2.73438, 3.125, 3.22266, 2.83203, 1.46484, 2.24609, 2.53906, 2.63672, 3.125, 3.51562, 2.63672, 1.17188, 2.63672, 2.53906, 2.53906, 2.92969, 3.125, 2.63672, 1.36719, 2.14844, 2.73438, 2.83203, 3.02734, 3.71094, 2.53906, 1.66016, 2.63672, 2.34375, 2.83203, 3.125, 3.02734, 2.63672, 1.5625, 2.34375, 2.83203, 2.53906, 3.22266, 3.41797, 2.63672, 1.36719, 2.14844, 2.53906, 2.53906, 3.125, 3.125, 2.63672, 1.5625, 2.14844, 2.63672, 2.83203, 3.125, 3.125, 2.92969, 1.26953, 2.14844, 2.44141, 2.63672, 3.41797, 3.125, 2.73438, 1.17188, 2.24609, 3.02734, 3.125, 3.41797, 2.73438, 1.85547, 2.73438, 3.02734, 2.63672, 2.53906, 4.6875, 2.73438, 3.22266, 1.07422, 0.878906, -0.0976562, -50.1953, 2.14844, 2.92969, 2.24609, 2.73438, 2.14844, 2.83203, 2.14844, 2.63672, 1.85547, 2.63672, 2.34375, 4.88281, 2.34375, 2.63672, 2.34375, 2.53906, 2.34375, 2.83203, 4.88281, 5.46875, 2.05078, 3.125, 2.53906};

  // Double_t lowDtCuts[numCombos] = {2.92969, 3.22266, 2.92969, 3.125,
  // 				   3.22266, 3.32031, 3.02734, 3.22266,
  // 				   3.41797, 3.51562, 2.92969, 3.125,
  // 				   3.02734, 3.22266, 3.41797, 3.22266};
  // Double_t highDtCuts[numCombos] = {5.46875, 5.85938, 5.37109, 5.85938,
  // 				    5.66406, 5.66406, 5.56641, 5.85938,
  // 				    5.66406, 5.76172, 5.56641, 5.76172,
  // 				    5.37109, 5.95703, 5.76172, 5.76172};

  // Double_t lowDtCuts[numCombos];
  // Double_t highDtCuts[numCombos];
  // for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
  //   lowDtCuts[comboInd] = -1000;
  //   highDtCuts[comboInd] = 1000;
  // }


  // Double_t lowDtCuts[numCombos] = {0};
  // Double_t highDtCuts[numCombos] = {0};  
  
  for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
    // Int_t ant1 = comboInd; // i.e. phi
    // // Int_t ant2 = (comboInd + 1)%NUM_PHI;
    // // Int_t ant2 = comboInd + 2*NUM_PHI;
    // Int_t ant2 = comboInd + NUM_PHI;
    // Int_t combo = cc->comboIndices[ant1][ant2];
    Int_t combo = comboInd;
    Int_t ant1 = cc->comboToAnt1s.at(comboInd);
    Int_t ant2 = cc->comboToAnt2s.at(comboInd);    
    // std:: cout << ant1 << "\t" << ant2 << "\t" << combo << std::endl;
    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
    
    TString name = TString::Format("hDtProf_%d_%d", ant1, ant2);
    TString title = TString::Format("2D profile of #deltat_{measured} as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);
    hDtProfs[comboInd] = new TProfile2D(name, title,
					numBinsPhi, phiDegMin, phiDegMax,
					numBinsTheta, thetaDegMin, thetaDegMax);

    name = TString::Format("hThetaPhiExpected_%d_%d", ant1, ant2);
    title = TString::Format("Number of events as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);

    const Int_t nDim1 = 2;
    Int_t nBins1[nDim1] = {numBinsPhi, numBinsTheta};
    Double_t mins1[nDim1] = {phiDegMin, thetaDegMin};
    Double_t maxs1[nDim1] = {phiDegMax, thetaDegMax};    
    
    hThetaPhiExpecteds[comboInd] = new THnSparseF(name,title, nDim1, nBins1, mins1, maxs1);

    name = TString::Format("hPhiExpPhiAnt_%d_%d", ant1, ant2);
    title = TString::Format("Antenna azimuth vs. expected azimuth for antennas %d and %d; Expected azimuth (degrees); Antenna azimuth (degrees)", ant1, ant2);
    const Int_t nDim2 = 2;
    Int_t nBins2[nDim2] = {numBinsPhi, numBinsTheta};
    Double_t mins2[nDim2] = {phiDegMin, phiDegMin};
    Double_t maxs2[nDim2] = {phiDegMax, phiDegMax};    

    hPhiExpPhiAnts[comboInd] = new THnSparseF(name, title, nDim2, nBins2, mins2, maxs2);


    name = TString::Format("hCorrDts_%d_%d", ant1, ant2); 
    title = TString::Format("Correlation coefficient vs #deltat_{measured} between antennas %d and %d; #deltat_{measured} (ns); Correlation coefficient (no units)", ant1, ant2);
    const Int_t nDim3 = 2;
    Int_t nBins3[nDim3] = {numDeltaTBins, numCorrBins};
    Double_t mins3[nDim3] = {minDeltaTBin, minCorrBin};
    Double_t maxs3[nDim3] = {maxDeltaTBin, maxCorrBin};    
    hCorrDts[comboInd] = new THnSparseF(name, title, nDim3, nBins3, mins3, maxs3);
    
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

	  const Int_t nDim = 2;
	  Double_t coords1[nDim] = {phiExpected, thetaExpected};
	  hThetaPhiExpecteds[comboInd]->Fill(coords1);

	  Double_t coords2[nDim] = {phiExpected, antPhiDeg1};
	  hPhiExpPhiAnts[comboInd]->Fill(coords2);

	  Double_t coords3[nDim] = {phiExpected, antPhiDeg2};	  
	  hPhiExpPhiAnts[comboInd]->Fill(coords3);	  
	}
	// }
	const Int_t nDim4 = 2;
	Double_t coords4[nDim4] = {correlationDeltaTs->at(pol).at(combo), correlationValues->at(pol).at(combo)};
	// std::cout << coords4[0] << "\t" << coords4[1] << std::endl;
	hCorrDts[comboInd]->Fill(coords4);
      }
    }
    p++;
  }


  std::vector<Double_t> minVals;
  std::vector<Double_t> maxVals;

  for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
    TH1D* hCorrDtProj = hCorrDts[comboInd]->Projection(0);
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

    minVals.push_back(minVal);
    maxVals.push_back(maxVal);    
    
    delete hCorrDtProj;
    
  }

  std::cout << "lowDtCuts[numCombos] = {";
  for(int comboInd=0; comboInd < numCombos; comboInd++){
    std::cout << minVals.at(comboInd);
    if(comboInd < numCombos-1){
      std::cout << ", ";
    }
  }
  std::cout << "};" << std::endl;

  std::cout << "highDtCuts[numCombos] = {";
  for(int comboInd=0; comboInd < numCombos; comboInd++){
    std::cout << maxVals.at(comboInd);
    if(comboInd < numCombos-1){
      std::cout << ", ";
    }
  }
  std::cout << "};" << std::endl;

  
  THnSparseF* hSparses[numCombos];
  //  const char* name, const char* title, Int_t dim, const Int_t* nbins, const Double_t* xmin = 0, const Double_t* xmax = 0, Int_t chunksize = 1024*16
  for(Int_t comboInd=0; comboInd<numCombos; comboInd++){
    Int_t ant1 = ant1s.at(comboInd);
    Int_t ant2 = ant2s.at(comboInd);

    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    TString title = TString::Format("2D profile of #deltat_{measured} as a function of expected azimuth and elevation for antennas %d and %d; Azimuth (degrees); Elevation (degrees)", ant1, ant2);

    const Int_t nDim = 2;
    Int_t nBins[nDim] = {numBinsPhi, numBinsTheta};
    Double_t mins[nDim] = {phiDegMin, thetaDegMin};
    Double_t maxs[nDim] = {phiDegMax, thetaDegMax};
    hSparses[comboInd] = new THnSparseF(name, title, nDim, nBins, mins, maxs);

    Double_t halfThetaBinWidth = 0.5*(hDtProfs[comboInd]->GetYaxis()->GetBinLowEdge(2) - hDtProfs[comboInd]->GetYaxis()->GetBinLowEdge(1));
    Double_t halfPhiBinWidth = 0.5*(hDtProfs[comboInd]->GetXaxis()->GetBinLowEdge(2) - hDtProfs[comboInd]->GetXaxis()->GetBinLowEdge(1));


    for(Int_t binx=1; binx<=hDtProfs[comboInd]->GetXaxis()->GetNbins(); binx++){
      Double_t phi = hDtProfs[comboInd]->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
      // phi *= TMath::DegToRad();
      for(Int_t biny=1; biny<=hDtProfs[comboInd]->GetYaxis()->GetNbins(); biny++){
	if(hDtProfs[comboInd]->GetBinContent(binx, biny) != 0){
	  Double_t theta = hDtProfs[comboInd]->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	  // theta *= TMath::DegToRad();
	  Double_t dt_m = hDtProfs[comboInd]->GetBinContent(binx, biny);

	  Double_t coords[nDim] = {phi, theta};
	  hSparses[comboInd]->Fill(coords, dt_m);
	}
      }
    }
    delete hDtProfs[comboInd];
    hSparses[comboInd]->Write();
    
    hThetaPhiExpecteds[comboInd]->Write();
    hPhiExpPhiAnts[comboInd]->Write();
    hCorrDts[comboInd]->Write();
    
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
