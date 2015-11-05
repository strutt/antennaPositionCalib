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



std::vector<THnSparseF*> hSparses;
std::vector<THnSparseF*> hSparses2;
AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArray;
std::vector<Double_t> phaseCentreToAmpaDeltaTs;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;

CrossCorrelator* cc;
const Int_t VARS_PER_ANT = 1; //3;
const Int_t numVars = NUM_SEAVEYS*VARS_PER_ANT;

Double_t sumOverSquaredDifferences(const Double_t* allTheVars);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTLookupHistogramsPlots.root");
  geom = AnitaGeomTool::Instance();
  geom->useKurtAnitaIIINumbers(1);
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
  // std::cout << geom->aftForeOffsetAngleVerticalKurtAnitaIII*TMath::RadToDeg() << std::endl;

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

  
  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   TString varName = TString::Format("deltaT_%d", ant);
  //   variables.at(varInd) = phaseCentreToAmpaDeltaTs.at(ant); //ns
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
  //   // if(ant==0){
  //   //   min->FixVariable(varInd);
  //   // }
  //   varInd++;
  // }
  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    TString varName = TString::Format("deltaR_%d", ant);
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    varInd++;
  }
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

  // for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
  //   phaseCentreToAmpaDeltaTs.at(ant) = vars[varInd];
  //   varInd++;
  //   // phaseCentreToAmpaDeltaTs.at(ant) = extraCableDelays[ant];    
  // }

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical]=rArray.at(ant)+vars[varInd];
    varInd++;
    // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant) + fittedDeltaRs[ant];    
  }
  
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

  for(UInt_t comboInd=0; comboInd<combos.size(); comboInd++){

    Int_t combo = combos.at(comboInd);
    THnSparseF* hSparse = hSparses.at(comboInd);
    THnSparseF* hSparse2 = hSparses2.at(comboInd);

    // std::cout << comboInd << "\t" << hSparse << "\t" << hSparse2 << std::endl;
    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);

    if(ant2 - ant1 == 16 || ant2 - ant1 == 32 || ((ant1/16) == (ant2/16) && (ant2 - ant1 == 1 || ant2 - ant1 == 15))){

    
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
	theta *= -1*TMath::DegToRad(); // INVERSION

	// This function doesn't care about offsets.
	// That only comes into play when considering a source location.
	Double_t dt_e = emptyPat.getDeltaTExpected(ant2, ant1, phi, theta);  

	// std::cout << dt_e << "\t" << dt_m << std::endl;

	// Is this the right way around?
	dt_m += phaseCentreToAmpaDeltaTs.at(ant1);
	dt_m -= phaseCentreToAmpaDeltaTs.at(ant2);

	if(dt_m_error > 0){
	  // Double_t diff = (dt_e - dt_m)/dt_m_error;
	  Double_t diff = (dt_e - dt_m);

	  // std::cout << diff << "\t" << (dt_e - dt_m) << "\t" << dt_m_error << std::endl;
	  
	  count++;
	  sumOfSquareDifferences += diff*diff;
	}
	  
	// if(phiBin==650 && thetaBin==87 && ant1==1 && ant2==2){
	//   std::cout << combo << "\t" << ant1 << "\t" << ant2 << "\t" << dt_e << "\t" 
	// 	    << dt_m << "\t" << diff << "\t" << geom->getAntR(ant1) << "\t"
	// 	    << geom->getAntR(ant2) << "\t" << phi << "\t" << theta << std::endl;
	// }      
      }
    }
  }
    // std::cout << TMath::Sqrt(sumOfSquareDifferences/count) << std::endl << std::endl;

  return TMath::Sqrt(sumOfSquareDifferences/count);
}


