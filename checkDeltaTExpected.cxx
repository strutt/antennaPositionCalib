// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to debug deltaT_expected in CrossCorrelator
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "TROOT.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"

int main(int argc, char *argv[])
{

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  TString lindaFileName = "newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, AnitaPol::kVertical);    
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }

  CrossCorrelator* cc = new CrossCorrelator();
  // cc->kDebug = 1;
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  gROOT->ProcessLine(".x /Users/benstrutt/Org/workNotes/weeklyMeeting/2016_02_29/meanDdtForAllAntennaPairs.C");
  TH2D* hMeanOffsets = (TH2D*) gDirectory->Get("h2");
  hMeanOffsets->SetName("hMeanOffsets");
    
  TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  lindaFileNameReference->Write();

  const Double_t minus6Deg = -6;

  Int_t ant1 = 1;
  Int_t ant2 = 2;
  
  // TH2D* hZoomed = cc->makeBlankZoomedImage("hTest", "Testing #deltat_{expected}",
  // 					   cc->phiArrayDeg[ant1].at(0), minus6Deg);
  
  TH2D* hZoomed = cc->makeZoomedImage(pol, 0xffff, cc->phiArrayDeg[ant1].at(0), minus6Deg);
  hZoomed->SetName("hTest");
  hZoomed->SetTitle("Testing #deltat_{expected}");
    

  for(int binx=1; binx<=hZoomed->GetNbinsX(); binx++){
    Double_t phiDeg = hZoomed->GetXaxis()->GetBinLowEdge(binx);
    Double_t phiWave = phiDeg*TMath::DegToRad();
    for(int biny=1; biny<=hZoomed->GetNbinsY(); biny++){
      Double_t thetaDeg = hZoomed->GetYaxis()->GetBinLowEdge(biny);
      Double_t thetaWave = thetaDeg*TMath::DegToRad();
      
      Double_t dt = cc->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);

      hZoomed->Fill(phiDeg+0.001, thetaDeg+0.001, dt);
    }
  }

  Double_t theMinus6DegBin = -1;
  Double_t deltaThetaDeg = 10000;
  for(int biny=1; biny<=hZoomed->GetNbinsY(); biny++){
    Double_t thetaDeg = hZoomed->GetYaxis()->GetBinLowEdge(biny);
    if(TMath::Abs(thetaDeg - minus6Deg) < deltaThetaDeg){
      deltaThetaDeg = TMath::Abs(thetaDeg - minus6Deg);
      theMinus6DegBin = biny;
    }
  }

  const Double_t DdtThresh = TMath::Abs(hMeanOffsets->GetBinContent(ant1+1, ant2+1));
  TH2D* hDeltaPhi = new TH2D("hDeltaPhi", "#delta#phi vs. #Delta#deltat; #delta#phi (Degrees); #Delta#deltat (ns);", 128, 0, 1, 128, 0, 2*DdtThresh);
  for(int binx=1; binx<=hZoomed->GetNbinsX(); binx++){
    Double_t phiDeg = hZoomed->GetXaxis()->GetBinLowEdge(binx);
    Double_t dt = hZoomed->GetBinContent(binx, theMinus6DegBin);
    for(int binx2=binx; binx2<=hZoomed->GetNbinsX(); binx2++){
      Double_t phiDeg2 = hZoomed->GetXaxis()->GetBinLowEdge(binx2);
      Double_t dt2 = hZoomed->GetBinContent(binx2, theMinus6DegBin);
      Double_t Ddt = TMath::Abs(dt2 - dt);
      
      if(Ddt >= DdtThresh){
	Double_t deltaPhiDeg = phiDeg2 - phiDeg;
	std::cout << binx << "\t" << phiDeg << "\t" << phiDeg2 << "\t" << deltaPhiDeg << "\t"
		  << dt << "\t" << dt2 << "\t" << Ddt << std::endl;

	hDeltaPhi->Fill(deltaPhiDeg, Ddt);
	break;
      }
    }
  }


  outFile->Write();
  outFile->Close();
  
  return 0;
}


