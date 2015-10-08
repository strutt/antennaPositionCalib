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

// const Int_t numCombos = NUM_PHI;
const Int_t numCombos = NUM_COMBOS;
THnSparseF* hSparses[numCombos];
TH2D* hProfs[numCombos];
AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArray;
std::vector<Double_t> phaseCentreToAmpaDeltaTs;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;

const Int_t VARS_PER_ANT = 4;
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
    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    hSparses[comboInd] = (THnSparseF*) fInputFile->Get(name);

    hProfs[comboInd] = NULL; 
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
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    if(ant==0){
      min->FixVariable(varInd);
    }
    varInd++;
    
    varName = TString::Format("deltaZ_%d", ant);
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    if(ant==0){
      min->FixVariable(varInd);
    }
    varInd++;

    varName = TString::Format("deltaPhi_%d", ant);
    variables.at(varInd) = 0;
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    if(ant==0){
      min->FixVariable(varInd);
    }
    varInd++;
    
    varName = TString::Format("phaseCentreToAmpaDeltaT_%d", ant);
    variables.at(varInd) = phaseCentreToAmpaDeltaTs.at(ant); //ns
    min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    if(ant==0){
      min->FixVariable(varInd);
    }
    varInd++;
  }

  

  assert(varInd==numVars);
  
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

    geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)+vars[varInd];
    varInd++;
    geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = zArray.at(ant)+vars[varInd];
    varInd++;
    geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = phiArray.at(ant)+vars[varInd];
    varInd++;
    phaseCentreToAmpaDeltaTs.at(ant) = vars[varInd];
    varInd++;
  }


  if(varInd != numVars){
    std::cerr << "Warning! Number of fitting variables mismatched!" << std::endl;
  }
  
  // for(int varInd=0; varInd<numVars; varInd++){
  //   Int_t ant = varInd + NUM_PHI;
  //   geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)+deltaR[varInd];
  // }
  
  for(int comboInd=0; comboInd<numCombos; comboInd++){

    // TH2D* prof = hProfs[comboInd];
    THnSparseF* hSparse = hSparses[comboInd];
    // std::cout << hSparse << std::endl;
    // Double_t halfThetaBinWidth = 0.5*(prof->GetYaxis()->GetBinLowEdge(2) - prof->GetYaxis()->GetBinLowEdge(1));
    // Double_t halfPhiBinWidth = 0.5*(prof->GetXaxis()->GetBinLowEdge(2) - prof->GetXaxis()->GetBinLowEdge(1));

    Double_t halfThetaBinWidth = 0.5*(hSparse->GetAxis(1)->GetBinLowEdge(2) - hSparse->GetAxis(1)->GetBinLowEdge(1));
    Double_t halfPhiBinWidth = 0.5*(hSparse->GetAxis(0)->GetBinLowEdge(2) - hSparse->GetAxis(0)->GetBinLowEdge(1));

    Int_t ant1 = ant1s.at(comboInd);
    Int_t ant2 = ant2s.at(comboInd);

    for(Int_t binInd=0; binInd < hSparse->GetNbins(); binInd++){
      Int_t coords[2] = {0};
      Double_t dt_m = hSparse->GetBinContent(binInd, coords);
      Int_t phiBin = coords[0];
      Int_t thetaBin = coords[1];
      Double_t theta = hSparse->GetAxis(1)->GetBinLowEdge(thetaBin) + halfThetaBinWidth;
      Double_t phi = hSparse->GetAxis(0)->GetBinLowEdge(phiBin) + halfPhiBinWidth;
      // std::cout << dt_m << '\t' << phiBin << "\t" << thetaBin << "\t" << phi << "\t" << theta << std::endl;
      
    //   std::cout << dt_m << "\t" << phiBin << "\t" << thetaBin << std::endl;
    // }
    
    
    // for(Int_t binx=1; binx<=prof->GetXaxis()->GetNbins(); binx++){
    //   Double_t phi = prof->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
      phi *= TMath::DegToRad();
      // for(Int_t biny=1; biny<=prof->GetYaxis()->GetNbins(); biny++){
      // 	if(prof->GetBinContent(binx, biny) != 0){
	  // Double_t theta = prof->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
      theta *= -1*TMath::DegToRad(); // INVERSION
      Double_t dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, theta);
      dt_e += phaseCentreToAmpaDeltaTs.at(ant1);
      dt_e -= phaseCentreToAmpaDeltaTs.at(ant2);
	  
      // Double_t dt_m = prof->GetBinContent(binx, biny);
      Double_t diff = dt_e - dt_m;
      count++;
      sumOfSquareDifferences += diff*diff;
      // 	}
      // }
    }
  }
  return TMath::Sqrt(sumOfSquareDifferences/count);
}



void makeHistosFromGeom(Int_t comboInd, TH1D* h1, TH2D* h2, TString nameAppended, TString titleAppended){


  if(hProfs[comboInd]==NULL){
    hProfs[comboInd] = hSparses[comboInd]->Projection(1, 0);
  }
  TH2D* prof = hProfs[comboInd];
  
  const Int_t numBinsPhi = prof->GetXaxis()->GetNbins();
  Double_t thetaDegMin = -50;
  Double_t thetaDegMax = 50;
  const Int_t numBinsTheta = prof->GetYaxis()->GetNbins();
  Double_t halfThetaBinWidth = 0.5*(prof->GetYaxis()->GetBinLowEdge(2) - prof->GetYaxis()->GetBinLowEdge(1));
  Double_t halfPhiBinWidth = 0.5*(prof->GetXaxis()->GetBinLowEdge(2) - prof->GetXaxis()->GetBinLowEdge(1));
  Double_t phiDegMin = 0;
  Double_t phiDegMax = 360;

  Int_t ant1 = ant1s.at(comboInd);
  Int_t ant2 = ant2s.at(comboInd);
  
  TString name = "hPureDiff_" + nameAppended + TString::Format("_%d_%d", ant1, ant2);
  TString title = TString::Format("Difference between measured and predicted #deltats for antennas %d and %d", ant1, ant2);
  title += titleAppended;
  title += "; #deltat_{expected} - #deltat_{measured} (ns)";
  title += "; Number of bins";

  if(h1!=NULL){
    delete h1;
  }
  if(h2!=NULL){
    delete h2;
  }
  
  h1 = new TH1D(name, title, 32, -0.5, 0.5);

  name = "hDtExpected_2D_" + nameAppended + TString::Format("_%d_%d", ant1, ant2);
  title = TString::Format("Predicted #deltats for antennas %d and %d", ant1, ant2);
  title += "; Azimuth (degrees)";
  title += "; Elevation (degrees)";  
  title += "; #delta t_{expected} (ns)";    
  
  h2 = new TH2D(name, title, 
		numBinsPhi, phiDegMin, phiDegMax,
		numBinsTheta, thetaDegMin, thetaDegMax);
    


  for(Int_t binx=1; binx<=prof->GetXaxis()->GetNbins(); binx++){
    Double_t phi = prof->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
    phi *= TMath::DegToRad();
    for(Int_t biny=1; biny<=prof->GetYaxis()->GetNbins(); biny++){
      if(prof->GetBinContent(binx, biny) != 0){
	Double_t theta = prof->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	theta *= TMath::DegToRad();

	Double_t dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, -theta);
	Double_t dt_m = prof->GetBinContent(binx, biny);

	//h1->Fill(dt_e); //-dt_m); //*(dt_e-dt_m));
	// std::cout << dt_e << "\t" << dt_m << std::endl;
	h1->Fill(dt_e-dt_m);

	// h2->SetBinContent(binx, biny, dt_e);
	h2->SetBinContent(binx, biny, dt_e-dt_m);	
      }
    }
  }
}
