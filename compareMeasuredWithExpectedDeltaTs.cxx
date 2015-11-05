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

#include <ProgressBar.h>
#include <CrossCorrelator.h>

const Double_t extraCableDelays[NUM_SEAVEYS] =
  {0,0.139733,0.00697198,0.0870043,0.0502864,0.139811,0.0564542,0.143131,
   0.0755872,0.10746, -0.0529778,0.109205,0.0199421,0.062739,0.0813765,
   0.033311,0.107176,0.141001,0.104494,0.16646,0.135572,0.14326,0.172411,
   0.113276,0.135898,0.0634696,0.078761,0.0890165,0.198665,0.0649681,0.057582,
   0.0800643,0.0635002,0.211805,0.133592,0.136322,0.125028,0.0841643,0.13517,
   0.0633426,0.0853349,0.0491331,0.0996015,0.0681319,0.0430019,0.0380842,
   0.0419707,-0.2944};

const Double_t fittedDeltaRs[NUM_SEAVEYS] =
  {-0.0959145, -0.0969593, -0.0976093, -0.0971427, -0.0972972, -0.0971141, -0.0971505, -0.0964771,
   -0.0960801, -0.094751 , -0.0938789, -0.0931123, -0.0941293, -0.0937031, -0.0944757, -0.0954,
   -0.0958407, -0.0967933, -0.0975012, -0.0970367, -0.0971045, -0.0971866, -0.0969005, -0.0965772,
   -0.0958794, -0.0949967, -0.0939146, -0.0936601, -0.0936972, -0.0940001, -0.0947252, -0.0954105,
   -0.0959116, -0.0966248, -0.0972106, -0.0971995, -0.0969414, -0.0970038, -0.0968214, -0.0965757,
   -0.0959988, -0.0949063, -0.0933371, -0.0934933, -0.0942268, -0.094005 , -0.0948142, -0.0957325};


std::vector<THnSparseF*> hSparses;
AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArray;
std::vector<Double_t> phaseCentreToAmpaDeltaTs;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;

const Int_t VARS_PER_ANT = 1;
const Int_t numVars = NUM_SEAVEYS*VARS_PER_ANT;

Double_t sumOverSquaredDifferences(const Double_t* allTheVars);
void makeHistosFromGeom(Int_t comboInd, TH1D* h1, TH2D* h2, TString nameAppended, TString titleAppended);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTLookupHistogramsPlots.root");
  CrossCorrelator* cc = new CrossCorrelator();

  for(Int_t comboInd=0; comboInd < NUM_COMBOS; comboInd++){
    // Int_t ant1 = comboInd; // i.e. phi
    // // Int_t ant2 = (comboInd + 1)%NUM_PHI;
    // // Int_t ant2 = comboInd + 2*NUM_PHI;
    // Int_t ant2 = comboInd + NUM_PHI;
    // Int_t combo = cc->comboIndices[ant1][ant2];

    Int_t combo = comboInd;
    Int_t ant1 = cc->comboToAnt1s.at(comboInd);
    Int_t ant2 = cc->comboToAnt2s.at(comboInd);

    if(!(TMath::Abs(ant1 - ant2) == NUM_PHI || TMath::Abs(ant1 - ant2) == (2*NUM_PHI) || TMath::Abs(ant1 - ant2) == 1 || TMath::Abs(ant1 - ant2) == (NUM_PHI-1))){
      continue;
    }
    
    // std:: cout << ant1 << "\t" << ant2 << "\t" << combo << std::endl;
    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);
    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    hSparses.push_back((THnSparseF*) fInputFile->Get(name));
  }
  
  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  
  geom = AnitaGeomTool::Instance();
  geom->fUseKurtAnitaIIINumbers = 1;
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArray.push_back(geom->getAntPhiPositionRelToAftFore(ant));
    phaseCentreToAmpaDeltaTs.push_back(0); //ns delay
    // std::cout << ant << "\t" << rArray.at(ant) << "\t" << zArray.at(ant)
    // 	      << "\t" << phiArray.at(ant)*TMath::RadToDeg() << std::endl;
  }

  // const Int_t numDefaults = 2;

  // TH2D* hDiff2D[numDefaults] = {NULL};
  // TH1D* hDiff[numDefaults] = {NULL};

  // geom->fUseKurtAnitaIIINumbers = 0;
  // makeHistosFromGeom(0, hDiff[0], hDiff2D[0], "feed", " for feed locations");
  // Double_t zeros1[1] = {0};
  // std::cout << sumOverSquaredDifferences(zeros1) << std::endl;
  // geom->fUseKurtAnitaIIINumbers = 1;
  // makeHistosFromGeom(0, hDiff[1], hDiff2D[1], "photo", " for photogrammetry numbers");
  // std::cout << sumOverSquaredDifferences(zeros1) << std::endl;
  
  // std::vector<Double_t> rSteps;
  // std::vector<Double_t> sumOverSquares;
  // Double_t myStepSize = 0.002;

  // Double_t minResidual = DBL_MAX;
  // Int_t bestRInd = 0;
  // const Int_t numDeltaRInds = 50;
  
  // for(int deltaRInd=0; deltaRInd < numDeltaRInds; deltaRInd++){
  //   // Double_t deltaR = 0.1*deltaRInd - 0.4;
  //   Double_t rStep = (myStepSize*deltaRInd - 0.2*myStepSize*numDeltaRInds);

  //   Double_t ddt2 = sumOverSquaredDifferences(&rStep);

  //   if(ddt2 < minResidual){
  //     bestRInd = deltaRInd;
  //     minResidual = ddt2;
  //   }

  //   // std::cout << geom << std::endl;
  //   // std::cout << rStep << "\t" << ddt2 << std::endl;
  //   rSteps.push_back(rStep);
  //   sumOverSquares.push_back(ddt2);
  // }

  // // makeHistosFromGeom(deltaRInd+2, hDiff[deltaRInd+2], hDiff2D[deltaRInd+2]);
  // // makeHistosFromGeom(bestRInd+2, hDiff[bestRInd+2], hDiff2D[bestRInd+2]);

  
  // TGraph* gr = new TGraph(rSteps.size(), &rSteps[0], &sumOverSquares[0]);
  // gr->SetName("grDeltaRScan");
  // gr->SetTitle("<#Deltat> for antennas 0 and 16; #deltar_{16} relative to photogrammetry positions (m); <#Deltat> (ns)");
  // gr->Write();






  




  
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
  // min->SetLimitedVariable(0, "deltaR", variables[0], step[0], 0, 10);

  Int_t varInd = 0;
  for(Int_t antInd=0; antInd < NUM_SEAVEYS; antInd++){
    Int_t ant = antInd;

    TString varName = TString::Format("deltaR_%d", ant);
    // variables.at(varInd) = 0;
    // min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    // // if(ant==0){
    // //   min->FixVariable(varInd);
    // // }
    // varInd++;
    
    varName = TString::Format("deltaZ_%d", ant);
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    if(ant==0){
      min->FixVariable(varInd);
    }
    varInd++;

    // varName = TString::Format("deltaPhi_%d", ant);
    // variables.at(varInd) = 0;
    // min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    // if(ant==0){
    //   min->FixVariable(varInd);
    // }
    // varInd++;
    
    // varName = TString::Format("phaseCentreToAmpaDeltaT_%d", ant);
    // variables.at(varInd) = phaseCentreToAmpaDeltaTs.at(ant); //ns
    // min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    // if(ant==0){
    //   min->FixVariable(varInd);
    // }
    // varInd++;
  }

  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }
  
  // for(Int_t varInd=0; varInd < numVars; varInd++){
  //   Int_t ant = NUM_PHI + varInd;
  //   TString varName = TString::Format("deltaR_%d", ant);
  //   variables.at(varInd) = rArray.at(ant);
  //   min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
  // }

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
  
  const double *xs = min->X();
  std::cout << "Minimum = " << min->MinValue() << std::endl;
  std::cout << "Minimum: f(" << xs[0] << ") = " << min->MinValue()  << std::endl;

  TGraph* grFitterMin = new TGraph();
  grFitterMin->SetPoint(0, xs[0], min->MinValue());
  grFitterMin->SetName("grFitterMin");
  grFitterMin->Write();
  
  outFile->Write();
  outFile->Close();

  delete cc;

  

  return 0;
}




Double_t sumOverSquaredDifferences(const Double_t* vars){

  Double_t sumOfSquareDifferences = 0;
  Int_t count = 0;

  Int_t varInd=0;

  for(Int_t antInd=0; antInd < NUM_SEAVEYS; antInd++){
    Int_t ant = antInd;

    // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)+vars[varInd];
    // varInd++;

    geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)+fittedDeltaRs[ant];
    
    geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = zArray.at(ant)+vars[varInd];
    varInd++;
    // geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = phiArray.at(ant)+vars[varInd];
    // varInd++;

    // phaseCentreToAmpaDeltaTs.at(ant) = vars[varInd];
    // varInd++;
    phaseCentreToAmpaDeltaTs.at(ant) = extraCableDelays[ant];    

  }


  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }
  
  // for(int varInd=0; varInd<numVars; varInd++){
  //   Int_t ant = varInd + NUM_PHI;
  //   geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)+deltaR[varInd];
  // }
  
  for(UInt_t comboInd=0; comboInd<combos.size(); comboInd++){

    THnSparseF* hSparse = hSparses.at(comboInd);
    Int_t ant1 = ant1s.at(comboInd);
    Int_t ant2 = ant2s.at(comboInd);

    for(Int_t binInd=0; binInd < hSparse->GetNbins(); binInd++){
      Int_t coords[2] = {0};
      Double_t dt_m = hSparse->GetBinContent(binInd, coords);
      Int_t phiBin = coords[0];
      Int_t thetaBin = coords[1];
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

      Double_t diff = dt_e - dt_m;
      count++;
      sumOfSquareDifferences += diff*diff;
    }
  }
  return TMath::Sqrt(sumOfSquareDifferences/count);
}
