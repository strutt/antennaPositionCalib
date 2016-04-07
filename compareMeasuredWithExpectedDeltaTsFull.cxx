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



AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArray;
std::vector<Double_t> phaseCentreToAmpaDeltaTs;
std::vector<Int_t> combos;

TTree* deltaTChain = NULL;
CrossCorrelator* cc;
const Int_t VARS_PER_ANT = 1; //3;
const Int_t numVars = NUM_SEAVEYS*VARS_PER_ANT;

UInt_t eventNumber = 0;
Double_t correlationDeltaTs[NUM_COMBOS] = {0};
// Double_t correlationValues[NUM_COMBOS] = {0};  
Double_t correlationDeltaTsClose[NUM_COMBOS] = {0};
// Double_t correlationValuesClose[NUM_COMBOS] = {0};  
Double_t thetaExpected = 0;
Double_t phiExpected = 0;
UsefulAdu5Pat usefulPat;


Double_t sumOverSquaredDifferences(const Double_t* allTheVars);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTTree_run352-352Plots.root");
  deltaTChain = (TTree*) fInputFile->Get("deltaTTree");

  deltaTChain->SetBranchAddress("eventNumber", &eventNumber);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTs[%d]", NUM_COMBOS), correlationDeltaTs);
  // deltaTChain->SetBranchAddress(TString::Format("correlationValues[%d]", NUM_COMBOS), correlationValuesClose);
  deltaTChain->SetBranchAddress(TString::Format("correlationDeltaTsClose[%d]", NUM_COMBOS), correlationDeltaTsClose);
  // deltaTChain->SetBranchAddress(TString::Format("correlationValuesClose[%d]", NUM_COMBOS), correlationValues);
  deltaTChain->SetBranchAddress("thetaExpected", &thetaExpected);
  deltaTChain->SetBranchAddress("phiExpected", &phiExpected);

  
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
  }
  
  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArray.push_back(geom->getAntPhiPosition(ant));
    if(phiArray.at(ant)<0){
      phiArray.at(ant)+=TMath::TwoPi();
    }
    if(phiArray.at(ant)>TMath::TwoPi()){
      phiArray.at(ant)-=TMath::TwoPi();
    }
    phaseCentreToAmpaDeltaTs.push_back(0); //ns delay
  }

  Double_t zeros[numVars] = {0};
  for(int varInd=0; varInd<numVars; varInd++){
    zeros[varInd] = 0;
  }
  std::cout << "Before fitting func = " << sumOverSquaredDifferences(zeros) << std::endl;
  
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

  
  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    TString varName = TString::Format("deltaT_%d", ant);
    variables.at(varInd) = phaseCentreToAmpaDeltaTs.at(ant); //ns
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    // if(ant==0){
    //   min->FixVariable(varInd);
    // }
    switch(ant){
      
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 18:
    case 19:
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    case 26:
    case 27:
    case 28:
    case 29:
    case 34:
    case 35:
    case 36:
    case 37:
    case 38:
    case 39:
    case 40:
    case 41:
    case 42:
    case 43:
    case 44:
    case 45:
      min->FixVariable(varInd);
      break;
    default:
      break;
    }
    
    varInd++;
  }
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   TString varName = TString::Format("deltaR_%d", ant);
  //   variables.at(varInd) = 0;
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
  //   switch(ant){
      
  //   case 2:
  //   case 3:
  //   case 4:
  //   case 5:
  //   case 6:
  //   case 7:
  //   case 8:
  //   case 9:
  //   case 10:
  //   case 11:
  //   case 12:
  //   case 13:
  //   case 18:
  //   case 19:
  //   case 20:
  //   case 21:
  //   case 22:
  //   case 23:
  //   case 24:
  //   case 25:
  //   case 26:
  //   case 27:
  //   case 28:
  //   case 29:
  //   case 34:
  //   case 35:
  //   case 36:
  //   case 37:
  //   case 38:
  //   case 39:
  //   case 40:
  //   case 41:
  //   case 42:
  //   case 43:
  //   case 44:
  //   case 45:
  //     min->FixVariable(varInd);
  //     break;
  //   default:
  //     break;
  //   }
  //   varInd++;
  // }
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   TString varName = TString::Format("deltaPhiDeg_%d", ant);
  //   variables.at(varInd) = 0;
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], TMath::RadToDeg()*step[varInd]);
  //   varInd++;
  // }
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   TString varName = TString::Format("deltaZ_%d", ant);
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
  
  outFile->Write();
  outFile->Close();

  delete cc;

  

  return 0;
}




Double_t sumOverSquaredDifferences(const Double_t* vars){

  Double_t sumOfSquareDifferences = 0;
  Int_t count = 0;


  Int_t varInd=0;

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    phaseCentreToAmpaDeltaTs.at(ant) = vars[varInd];
    varInd++;
    // phaseCentreToAmpaDeltaTs.at(ant) = extraCableDelays[ant];    
  }

  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical]=rArray.at(ant)+vars[varInd];
  //   varInd++;
  //   // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant) + fittedDeltaRs[ant];    
  // }
  
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = phiArray.at(ant)+(TMath::DegToRad()*vars[varInd]);
  //   varInd++;
  //   // geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = phiArray.at(ant)+TMath::DegToRad()*fittedDeltaPhiDeg[ant];
  // }

  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){  
  //   geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = zArray.at(ant)+vars[varInd];
  //   varInd++;
  //   // geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = zArray.at(ant)+fittedDeltaZ[ant];
  // }

  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }

  // for(Long64_t entry=0; entry<deltaTChain->GetEntries(); entry++){
  deltaTChain->GetEntry(0);

    // if(eventNumber != 60832108) continue;
    std::cout << eventNumber << std::endl;
    for(UInt_t comboInd=0; comboInd<combos.size(); comboInd++){

      Int_t combo = combos.at(comboInd);

      // std::cout << comboInd << "\t" << hSparse << "\t" << hSparse2 << std::endl;
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      if(correlationDeltaTs[combo]==correlationDeltaTsClose[combo]){
      
	Double_t dt_e = usefulPat.getDeltaTExpected(ant2, ant1,
						    phiExpected*TMath::DegToRad(),
						    thetaExpected*TMath::DegToRad());
	Double_t dt_m = correlationDeltaTs[combo];

	// Is this the right way around?
	dt_m += phaseCentreToAmpaDeltaTs.at(ant1);
	dt_m -= phaseCentreToAmpaDeltaTs.at(ant2);
	
	Double_t diff = (dt_m - dt_e);
	sumOfSquareDifferences += diff*diff;
	count++;
      }
    }
  // }
  std::cout << TMath::Sqrt(sumOfSquareDifferences/count) << std::endl;
  return TMath::Sqrt(sumOfSquareDifferences/count);
}


