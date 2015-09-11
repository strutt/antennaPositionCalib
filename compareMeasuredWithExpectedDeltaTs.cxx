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

TProfile2D* prof;
AnitaGeomTool* geom;
Int_t ant1 = 0;
Int_t ant2 = 16;
UsefulAdu5Pat emptyPat;
Double_t halfThetaBinWidth;
Double_t halfPhiBinWidth;
CrossCorrelator* cc;

Double_t sumOverSquaredDifferences(const Double_t* allTheVars);
void makeHistosFromGeom(Int_t geomIndex, TH1D*& h1, TH2D*& h2);

int main(int argc, char *argv[])
{

  if(argc==0){
    std::cerr << "This can't happen" << std::endl;
  }
  
  TFile* fInputFile = TFile::Open("generateDeltaTLookupHistogramsPlots.root");
  prof = (TProfile2D*) fInputFile->Get("hProf");

  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");
  
  // OK sweet, now let's get deltaTExpected for the Prioritizerd and the photogrammetry positions  
  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  halfThetaBinWidth = 0.5*(prof->GetYaxis()->GetBinLowEdge(2) - prof->GetYaxis()->GetBinLowEdge(1));
  halfPhiBinWidth = 0.5*(prof->GetXaxis()->GetBinLowEdge(2) - prof->GetXaxis()->GetBinLowEdge(1));


  geom = AnitaGeomTool::Instance();
  // std::cout << geom->getAntR(ant1) << "\t"
  // 	    << geom->getAntPhiPositionRelToAftFore(ant1) << "\t"
  // 	    << geom->getAntZ(ant1) << std::endl;
  // std::cout << geom->getAntR(ant2) << "\t"
  // 	    << geom->getAntPhiPositionRelToAftFore(ant2) << "\t"
  // 	    << geom->getAntZ(ant2) << std::endl;
  // std::cout << std::endl << std::endl;

  cc = new CrossCorrelator();

  const Int_t numDefaults = 3;
  
  TH2D* hDiff2D[numDefaults] = {NULL};
  TH1D* hDiff[numDefaults] = {NULL};  


  
  
  
  for(Int_t defaultPosition=0; defaultPosition < numDefaults-1; defaultPosition++){
    makeHistosFromGeom(defaultPosition, hDiff[defaultPosition], hDiff2D[defaultPosition]);
  }
  
  delete cc;

  // Right, now let's try to minimize this bastard...
  // Mostly copied from the ROOT Minuit tutorial 
  const Int_t numVars = 6;
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor funcToMin(&sumOverSquaredDifferences, numVars);

  Double_t stepSize = 1e-3;
  std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

  // starting point
  std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);

  variables.at(0) = geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical];
  variables.at(1) = geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];
  variables.at(2) = geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical];
  variables.at(3) = geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];
  variables.at(4) = geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical];
  variables.at(5) = geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical];

  min->SetFunction(funcToMin);

  // Set the free variables to be minimized!
  min->SetVariable(0, "r1", variables[0], step[0]);
  min->SetVariable(1, "r2", variables[1], step[1]);
  min->SetVariable(2, "z1", variables[2], step[2]);
  min->SetVariable(3, "z2", variables[3], step[3]);
  min->SetVariable(4, "phi1", variables[4], step[4]);
  min->SetVariable(5, "phi2", variables[5], step[5]);

  min->FixVariable(0);
  min->FixVariable(2);
  min->FixVariable(4);  
  
  // do the minimization
  min->Minimize(); 
 
  const double *xs = min->X();
  std::cout << "Minimum = " << min->MinValue() << std::endl;

  geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[0];
  geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[1];
  geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[2];
  geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[3];
  geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[4];
  geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = xs[5];

  geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[0];
  geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[2];
  geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = xs[4];
  
  makeHistosFromGeom(2, hDiff[2], hDiff2D[2]);
  
  outFile->Write();
  outFile->Close();

  return 0;
}




Double_t sumOverSquaredDifferences(const Double_t* allTheVars){

  geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = allTheVars[0];
  geom->rPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = allTheVars[1];
  geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = allTheVars[2];
  geom->zPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = allTheVars[3];
  geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant1][AnitaPol::kVertical] = allTheVars[4];
  geom->azPhaseCentreFromVerticalHornKurtAnitaIII[ant2][AnitaPol::kVertical] = allTheVars[5];
  
  Double_t sumOfSquareDifferences = 0;
  for(Int_t binx=1; binx<=prof->GetXaxis()->GetNbins(); binx++){
    Double_t phi = prof->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
    phi *= TMath::DegToRad();
    for(Int_t biny=1; biny<=prof->GetYaxis()->GetNbins(); biny++){
      if(prof->GetBinContent(binx, biny) > 0){
	Double_t theta = prof->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	theta *= TMath::DegToRad();
	Double_t dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, -theta);
	Double_t dt_m = prof->GetBinContent(binx, biny);
	sumOfSquareDifferences += (dt_e-dt_m)*(dt_e-dt_m);
      }
    }
  }
  return sumOfSquareDifferences;
}



void makeHistosFromGeom(Int_t defaultPosition, TH1D*& h1, TH2D*& h2){

  const Int_t numBinsPhi = prof->GetXaxis()->GetNbins();
  Double_t thetaDegMin = -7;
  Double_t thetaDegMax = -5;
  const Int_t numBinsTheta = prof->GetYaxis()->GetNbins();
  halfThetaBinWidth = 0.5*(prof->GetYaxis()->GetBinLowEdge(2) - prof->GetYaxis()->GetBinLowEdge(1));
  halfPhiBinWidth = 0.5*(prof->GetXaxis()->GetBinLowEdge(2) - prof->GetXaxis()->GetBinLowEdge(1));
  Double_t phiDegMin = 0;
  Double_t phiDegMax = 360;

  
  TString name = TString::Format("hPureDiff%d", defaultPosition);
  TString title = "Difference between measured and predicted #deltats";
  switch (defaultPosition){
  case 0:
    title += " for Feed locations";
    break;
  case 1:
    title += " for photogrammetry locations ";
    break;
  case 2:
    title += " for optimized locations ";
    break;
    
  default:
    std::cerr << "Uh oh!" << std::endl;
    break;
  }

  h1 = new TH1D(name, title, 32, -0.5, 0.5);

  name += "2D";
  h2 = new TH2D(name, title, 
		numBinsPhi, phiDegMin, phiDegMax,
		numBinsTheta, thetaDegMin, thetaDegMax);
    
  for(Int_t binx=1; binx<=prof->GetXaxis()->GetNbins(); binx++){
    Double_t phi = prof->GetXaxis()->GetBinLowEdge(binx) + halfPhiBinWidth;
    phi *= TMath::DegToRad();
    for(Int_t biny=1; biny<=prof->GetYaxis()->GetNbins(); biny++){
      if(prof->GetBinContent(binx, biny) > 0){
	Double_t theta = prof->GetYaxis()->GetBinLowEdge(biny) + halfThetaBinWidth;
	theta *= TMath::DegToRad();

	Double_t dt_e = 0;
	switch (defaultPosition){
	case 0:
	  dt_e = cc->getDeltaTExpectedDouble(ant2, ant1, phi, theta);
	  break;	    
	default:
	  dt_e = emptyPat.getDeltaTExpected(ant1, ant2, phi, -theta);
	  break;
	}

	Double_t dt_m = prof->GetBinContent(binx, biny);

	h1->Fill(dt_e-dt_m); //*(dt_e-dt_m));
	h2->SetBinContent(binx, biny, (dt_e-dt_m)); //*(dt_e-dt_m));	  
      }
    }
  }
}
