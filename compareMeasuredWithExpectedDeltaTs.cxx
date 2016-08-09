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
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <THnSparse.h>

#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>
#include "AnitaEventCalibrator.h"

#include <ProgressBar.h>
#include <CrossCorrelator.h>
#include "OutputConvention.h"

TString variableNamePrefix;
std::vector<THnSparseF*> hSparses;
std::vector<THnSparseF*> hSparses2;
AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArrayDeg;
std::vector<Double_t> phaseCentreToAmpaDeltaTs;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;
Bool_t fMakePlots;

CrossCorrelator* cc;
const Int_t VARS_PER_ANT = 1; //3;
const Int_t numVars = NUM_SEAVEYS*VARS_PER_ANT;
const AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
Double_t extraCableDelays[NUM_SEAVEYS] = {-0.0610646, 0.0814919, -0.0318929, 0.0213868, 0.00642687, 0.0739681, -0.0187586, 0.038752, -0.0471173, -0.0213559, -0.196427, -0.0150133, -0.0971924, -0.0465036, -0.0127587, -0.0574128, 0.0104067, 0.0426622, 0.00792418, 0.0659761, 0.0458048, 0.0615984, 0.0971525, 0.0355208, 0.0565163, -0.0156482, 0.00579238, 0.0171781, 0.120797, -0.0179532, -0.0289847, -0.0103913, -0.0442358, 0.104188, 0.0176335, 0.0288141, 0.0228393, -0.00693136, 0.0572457, -0.0108378, 0.0149809, -0.0182012, 0.030143, 0.00987049, -0.0159015, -0.0346651, -0.0411542, -0.388232};

Double_t sumOverSquaredDifferences(const Double_t* allTheVars);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTLookupHistogramsPlots.root");
  geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(int surf=0; surf<12; surf++){
    for(int chan=0; chan < 9; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
    }
  }
  
  cc = new CrossCorrelator();

  for(Int_t comboInd=0; comboInd < NUM_COMBOS; comboInd++){
    // Int_t ant1 = comboInd; // i.e. phi
    // // Int_t ant2 = (comboInd + 1)%NUM_PHI;
    // // Int_t ant2 = comboInd + 2*NUM_PHI;
    // Int_t ant2 = comboInd + NUM_PHI;
    // Int_t combo = cc->comboIndices[ant1][ant2];

    Int_t combo = comboInd;
    Int_t ant1 = cc->comboToAnt1s.at(comboInd);
    Int_t ant2 = cc->comboToAnt2s.at(comboInd);

    if(!(TMath::Abs(ant1 - ant2) == NUM_PHI || TMath::Abs(ant1 - ant2) == (2*NUM_PHI) || (ant1/NUM_PHI==ant2/NUM_PHI && (TMath::Abs(ant1 - ant2) == 1 || TMath::Abs(ant1 - ant2) == (NUM_PHI-1))))){
      continue;
    }
    
    // std:: cout << ant1 << "\t" << ant2 << "\t" << combo << std::endl;
    combos.push_back(combo);
    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    hSparses.push_back((THnSparseF*) fInputFile->Get(name));

    name = TString::Format("hDtErrorSparse_%d_%d", ant1, ant2);
    hSparses2.push_back((THnSparseF*) fInputFile->Get(name));
  }
  
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArrayDeg.push_back(geom->getAntPhiPositionRelToAftFore(ant)*TMath::RadToDeg());
    phaseCentreToAmpaDeltaTs.push_back(0); //ns delay
  }
  // std::cout << geom->aftForeOffsetAngleVerticalKurtAnitaIII*TMath::RadToDeg() << std::endl;

  Double_t zeros[numVars] = {0};
  for(int varInd=0; varInd<numVars; varInd++){
    zeros[varInd] = 0;
  }
  
  // Right, now let's try to minimize this bastard...
  // Mostly copied from the ROOT Minuit tutorial 

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
   min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000);  // for GSL
   min->SetTolerance(0.0001);
   min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor funcToMin(&sumOverSquaredDifferences, numVars);

  Double_t stepSize = 1e-3;
  std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

  // starting point
  std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);

  min->SetFunction(funcToMin);

  // Set the free variables to be minimized!

  // const Double_t minDeltaPhiDeg = -5;
  // const Double_t maxDeltaPhiDeg = 5;
  
  Int_t varInd = 0;
  
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   variableNamePrefix = "deltaT";
  //   TString varName = variableNamePrefix + TString::Format("_%d", ant);
  //   variables.at(varInd) = phaseCentreToAmpaDeltaTs.at(ant); //ns
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
  //   varInd++;
  // }
  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    variableNamePrefix = "deltaR";
    TString varName = variableNamePrefix + TString::Format("_%d", ant);
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    varInd++;
  }
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   variableNamePrefix = "deltaPhiDeg";
  //   TString varName = variableNamePrefix + TString::Format("_%d", ant);
  //   variables.at(varInd) = 0;
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], TMath::RadToDeg()*step[varInd]);
  //   varInd++;
  // }
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   variableNamePrefix = "deltaZ";
  //   TString varName = variableNamePrefix + TString::Format("_%d", ant);
  //   variables.at(varInd) = 0;
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
  //   if(ant==0){
  //     min->FixVariable(varInd);
  //   }
  //   varInd++;
  // }

  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }

  TString outFileName = TString::Format("%s_", argv[0]) + variableNamePrefix + "Plots.root";
  TFile* outFile = new TFile(outFileName, "recreate");
  fMakePlots = true;
  TString variableNamePrefixTemp = variableNamePrefix;
  variableNamePrefix = "before";
  std::cout << "Before fitting func = " << sumOverSquaredDifferences(zeros) << std::endl;
  variableNamePrefix = variableNamePrefixTemp;
  fMakePlots = false;

  // Time it
  TStopwatch watch;
  watch.Start(kTRUE);

  // do the minimization
  min->Minimize(); 

  // Time!
  watch.Start(kFALSE);
  Int_t seconds = Int_t(watch.RealTime());
  Int_t hours = seconds / 3600;
  hours = hours < 0 ? 0 : hours;
  seconds = seconds - hours * 3600;
  Int_t mins = seconds / 60;
  mins = mins < 0 ? 0 : mins;
  seconds = seconds - mins * 60;
  fprintf(stderr, "Minimization took %02d:%02d:%02d\n", hours, mins, seconds);
  
  std::cout << "Minimum = " << min->MinValue() << std::endl;

  fMakePlots = true;
  // std::vector<Double_t> fittedVariables(numVars);
  // for(Int_t varInd2=0; varInd2<numVars; varInd2++){
  //   fittedVariables.at(varInd2) = min->X()[varInd2];
  // }
  sumOverSquaredDifferences(min->X());
  
  outFile->Write();
  outFile->Close();

  delete cc;

  

  return 0;
}




Double_t sumOverSquaredDifferences(const Double_t* vars){

  Double_t sumOfSquareDifferences = 0;
  Int_t comboCount = 0;

  Int_t varInd=0;

  TGraph* grs[NUM_COMBOS] = {NULL};

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    // phaseCentreToAmpaDeltaTs.at(ant) = vars[varInd];
    // varInd++;
    phaseCentreToAmpaDeltaTs.at(ant) = extraCableDelays[ant];    
  }

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    // cc->rArray.at(ant) = rArray.at(ant)+vars[varInd];
    cc->rArray[pol].at(ant) = rArray.at(ant)+vars[varInd];    
    varInd++;
  }
 
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   cc->phiArrayDeg.at(ant) = phiArrayDeg.at(ant) + vars[varInd];
  //   varInd++;
  // }

  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   cc->zArray.at(ant) = zArray.at(ant) + vars[varInd];
  //   varInd++;
  // }

  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }

  // Int_t comboCounter2[NUM_SEAVEYS] = {0};
  for(UInt_t comboInd=0; comboInd<combos.size(); comboInd++){

    Int_t combo = combos.at(comboInd);
    THnSparseF* hSparse = hSparses.at(comboInd);
    THnSparseF* hSparse2 = hSparses2.at(comboInd);

    // std::cout << comboInd << "\t" << hSparse << "\t" << hSparse2 << std::endl;
    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);

    if(ant2 - ant1 == 16 || ant2 - ant1 == 32 || ((ant1/16) == (ant2/16) && (ant2 - ant1 == 1 || ant2 - ant1 == 15))){

      if(fMakePlots){
	grs[combo] = new TGraph();
      }
      TGraph* gr = grs[combo];
      
      // comboCounter2[ant1]++;
      // comboCounter2[ant2]++;
      
      Double_t littleSumSquare = 0;
      Double_t binCount = 0;  
      for(Int_t binInd=0; binInd < hSparse->GetNbins(); binInd++){
	Int_t coords[2] = {0};
	Int_t coords2[2] = {0};	
	
	Double_t dt_m = hSparse->GetBinContent(binInd, coords);

	Int_t phiBin = coords[0];
	Int_t thetaBin = coords[1];

	Double_t dt_m_error = hSparse2->GetBinContent(binInd, coords2);

	// if(dt_m_error <= 0){
	//   std::cout << ant1 << "\t" << ant2 << "\t" << dt_m_error << "\t" << dt_m << std::endl;
	//   std::cout << coords[0] << "\t" << coords2[0] << "\t"
	// 	    << coords[1] << "\t" << coords2[1] << std::endl;
	// }
	
	// But these perhaps have been affected by rotation
	Double_t theta = hSparse->GetAxis(1)->GetBinLowEdge(thetaBin);
	Double_t phi = hSparse->GetAxis(0)->GetBinLowEdge(phiBin);

	phi *= TMath::DegToRad();
	// theta *= -1*TMath::DegToRad(); // INVERSION
	theta *= TMath::DegToRad();

	// This function doesn't care about offsets.
	// That only comes into play when considering a source location.
	// Double_t dt_e = emptyPat.getDeltaTExpected(ant2, ant1, phi, theta);
	Double_t dt_e = cc->getDeltaTExpected(pol, ant1, ant2, phi, theta);

	// std::cout << dt_e << "\t" << dt_m << std::endl;

	// Is this the right way around?
	dt_m += phaseCentreToAmpaDeltaTs.at(ant1);
	dt_m -= phaseCentreToAmpaDeltaTs.at(ant2);

	if(dt_m_error > 0){
	  Double_t diff = (dt_e - dt_m);
	  if(gr){
	    gr->SetPoint(gr->GetN(), phi*TMath::RadToDeg(), diff);
	  }
	  
	  // std::cout << diff << "\t" << (dt_e - dt_m) << "\t" << dt_m_error << std::endl;
	  littleSumSquare += 1e12*diff*diff;
	  binCount++;
	}
      }
      littleSumSquare/=binCount;
      // std::cout << combo << "\t" << ant1 << "\t" << ant2 << "\t" << binCount << "\t" << littleSumSquare << std::endl;
      sumOfSquareDifferences += littleSumSquare;
      comboCount++;

      if(fMakePlots){
	gr->SetName(TString::Format("gr_%d_%s", combo, variableNamePrefix.Data()));
	gr->SetTitle(TString::Format("ant1 = %d, ant2 = %d", ant1, ant2));
	gr->Write();
      }
      
    }
  }

  
  return TMath::Sqrt(sumOfSquareDifferences/comboCount);
}


