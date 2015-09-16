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


#include <RawAnitaHeader.h>
#include <UsefulAdu5Pat.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

#include <ProgressBar.h>
#include <CrossCorrelator.h>

const Int_t numCombos = 1; //NUM_PHI;
TProfile2D* profs[numCombos];
AnitaGeomTool* geom;
UsefulAdu5Pat emptyPat;
std::vector<Double_t> rArray;
std::vector<Double_t> zArray;
std::vector<Double_t> phiArray;
std::vector<Int_t> combos;
std::vector<Int_t> ant1s;
std::vector<Int_t> ant2s;

Double_t sumOverSquaredDifferences(const Double_t* allTheVars);
void makeHistosFromGeom(Int_t geomIndex, TH1D*& h1, TH2D*& h2);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTLookupHistogramsPlots.root");
  CrossCorrelator* cc = new CrossCorrelator();

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
    profs[comboInd] = (TProfile2D*) fInputFile->Get(name);
  }
  

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  


  geom = AnitaGeomTool::Instance();
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArray.push_back(geom->getAntPhiPositionRelToAftFore(ant));
  }

  const Int_t numDeltaRInds = 40;
  // const Int_t numDefaults = 3 + numDeltaRInds;

  // TH2D* hDiff2D[numDefaults] = {NULL};
  // TH1D* hDiff[numDefaults] = {NULL};

  // geom->fUseKurtAnitaIIINumbers = 0;
  // makeHistosFromGeom(0, hDiff[0], hDiff2D[0]);
  // geom->fUseKurtAnitaIIINumbers = 1;
  // makeHistosFromGeom(1, hDiff[1], hDiff2D[1]);

  std::vector<Double_t> facts;
  std::vector<Double_t> sumOverSquares;  
  Double_t myStepSize = 0.002;

  Double_t minResidual = DBL_MAX;
  Int_t bestRInd = 0;
  
  for(int deltaRInd=0; deltaRInd < numDeltaRInds; deltaRInd++){
    // Double_t deltaR = 0.1*deltaRInd - 0.4;
    Double_t fact = 1 + (myStepSize*deltaRInd - 0.5*myStepSize*numDeltaRInds);

    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant][AnitaPol::kVertical] = rArray.at(ant)*fact;
    }
    
    Double_t zeros[1] = {0};
    Double_t ddt2 = sumOverSquaredDifferences(zeros);

    if(ddt2 < minResidual){
      bestRInd = deltaRInd;
      minResidual = ddt2;
    }

    // std::cout << geom << std::endl;
    std::cout << fact << "\t" << ddt2 << std::endl;
    facts.push_back(fact);
    sumOverSquares.push_back(ddt2);
  }

  // makeHistosFromGeom(deltaRInd+2, hDiff[deltaRInd+2], hDiff2D[deltaRInd+2]);
  // makeHistosFromGeom(bestRInd+2, hDiff[bestRInd+2], hDiff2D[bestRInd+2]);  
  
  TGraph* gr = new TGraph(facts.size(), &facts[0], &sumOverSquares[0]);
  gr->SetName("grFacts");
  gr->SetTitle("Mean (#deltat_{measured} - #deltat_{expected})^{2} for antennas 0 and 16; Radial scale factor relative to photogrammetry positions; Mean (#deltat_{measured} - #deltat_{expected})^{2} (ns^{2})");
  gr->Write();
  
  // // Right, now let's try to minimize this bastard...
  // // Mostly copied from the ROOT Minuit tutorial 
  // const Int_t numVars = 1;
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // // set tolerance , etc...
  // min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  // min->SetMaxIterations(10000);  // for GSL 
  // min->SetTolerance(0.0001);
  // min->SetPrintLevel(1);

  // // create funciton wrapper for minmizer
  // // a IMultiGenFunction type 
  // ROOT::Math::Functor funcToMin(&sumOverSquaredDifferences, numVars);

  // Double_t stepSize = 1e-3;
  // std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

  // // starting point
  // std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);

  // variables.at(0) = 0; 
  // // variables.at(1) = geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];
  // // variables.at(2) = geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical];
  // // variables.at(3) = geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];
  // // variables.at(4) = geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical];
  // // variables.at(5) = geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];

  // min->SetFunction(funcToMin);

  // // Set the free variables to be minimized!
  // // min->SetLimitedVariable(0, "deltaR", variables[0], step[0], 0, 10);
  // min->SetVariable(0, "deltaR", variables[0], step[0]);
  // // min->SetVariable(1, "r2", variables[1], step[1]);
  // // min->SetVariable(2, "z1", variables[2], step[2]);
  // // min->SetVariable(3, "z2", variables[3], step[3]);
  // // min->SetVariable(4, "phi1", variables[4], step[4]);
  // // min->SetVariable(5, "phi2", variables[5], step[5]);

  // // min->FixVariable(2);
  // // min->FixVariable(3);  
  // // min->FixVariable(4);
  // // min->FixVariable(5);  
  
  // // do the minimization
  // min->Minimize(); 
 
  // const double *xs = min->X();
  // std::cout << "Minimum = " << min->MinValue() << std::endl;

  // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = rArray.at(ant1) + xs[0];
  // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = rArray.at(ant2) + xs[0];
  
  // // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[0];
  // // geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[2];
  // // geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[4];

  // // geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[1];
  // // geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[3];
  // // geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[5];

  // makeHistosFromGeom(numDefaults-1, hDiff[numDefaults-1], hDiff2D[numDefaults-1]);

  // TProfile2D* hpMeas = (TProfile2D*) prof->Clone("hpMeas");
  // hpMeas->SetTitle("#deltat_{measured}");  
  outFile->Write();
  outFile->Close();

  delete cc;

  

  return 0;
}




Double_t sumOverSquaredDifferences(const Double_t* allTheVars){

  Double_t sumOfSquareDifferences = allTheVars[0] + 0;
  Int_t count = 0;
  
  for(int comboInd=0; comboInd<numCombos; comboInd++){
    Double_t halfThetaBinWidth = 0.5*(profs[comboInd]->GetYaxis()->GetBinLowEdge(2) - profs[comboInd]->GetYaxis()->GetBinLowEdge(1));
    Double_t halfPhiBinWidth = 0.5*(profs[comboInd]->GetXaxis()->GetBinLowEdge(2) - profs[comboInd]->GetXaxis()->GetBinLowEdge(1));

    Int_t ant1 = ant1s.at(comboInd);
    Int_t ant2 = ant2s.at(comboInd);    
    
    for(Int_t binx=1; binx<=profs[comboInd]->GetXaxis()->GetNbins(); binx++){
      Double_t phi = profs[comboInd]->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
      phi *= TMath::DegToRad();
      for(Int_t biny=1; biny<=profs[comboInd]->GetYaxis()->GetNbins(); biny++){
	if(profs[comboInd]->GetBinContent(binx, biny) != 0){
	  Double_t theta = profs[comboInd]->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	  theta *= -1*TMath::DegToRad(); // INVERSION
	  Double_t dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, theta);
	  Double_t dt_m = profs[comboInd]->GetBinContent(binx, biny);
	  Double_t diff = dt_e - dt_m;
	  count++;
	  sumOfSquareDifferences += diff*diff;
	}
      }
    }
  }
  return sumOfSquareDifferences/count;
}



// void makeHistosFromGeom(Int_t defaultPosition, TH1D*& h1, TH2D*& h2){

//   const Int_t numBinsPhi = prof->GetXaxis()->GetNbins();
//   Double_t thetaDegMin = -50;
//   Double_t thetaDegMax = 50;
//   const Int_t numBinsTheta = prof->GetYaxis()->GetNbins();
//   halfThetaBinWidth = 0.5*(prof->GetYaxis()->GetBinLowEdge(2) - prof->GetYaxis()->GetBinLowEdge(1));
//   halfPhiBinWidth = 0.5*(prof->GetXaxis()->GetBinLowEdge(2) - prof->GetXaxis()->GetBinLowEdge(1));
//   Double_t phiDegMin = 0;
//   Double_t phiDegMax = 360;

  
//   TString name = TString::Format("hPureDiff%d", defaultPosition);
//   TString title = "Difference between measured and predicted #deltats";
//   switch (defaultPosition){
//   case 0:
//     title += " for Feed locations";
//     break;
//   case 1:
//     title += " for photogrammetry locations ";
//     break;
//   case 2:
//     title += " for optimized locations ";
//     break;
    
//   default:
//     // std::cerr << "Uh oh!" << std::endl;
//     break;
//   }

//   h1 = new TH1D(name, title, 32, -0.5, 0.5);

//   name += "_2D";
//   h2 = new TH2D(name, title, 
// 		numBinsPhi, phiDegMin, phiDegMax,
// 		numBinsTheta, thetaDegMin, thetaDegMax);
    
//   for(Int_t binx=1; binx<=prof->GetXaxis()->GetNbins(); binx++){
//     Double_t phi = prof->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
//     phi *= TMath::DegToRad();
//     for(Int_t biny=1; biny<=prof->GetYaxis()->GetNbins(); biny++){
//       if(prof->GetBinContent(binx, biny) != 0){
// 	Double_t theta = prof->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
// 	theta *= TMath::DegToRad();

// 	Double_t dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, -theta);
// 	Double_t dt_m = prof->GetBinContent(binx, biny);

// 	//h1->Fill(dt_e); //-dt_m); //*(dt_e-dt_m));
// 	// std::cout << dt_e << "\t" << dt_m << std::endl;
// 	h1->Fill(dt_e-dt_m);

// 	// h2->SetBinContent(binx, biny, dt_e);
// 	h2->SetBinContent(binx, biny, dt_e-dt_m);	
//       }
//     }
//   }
// }
