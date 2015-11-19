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


UsefulAdu5Pat usefulPat;
TProfile* hDeltaThetaDeg_pfx;
TProfile* hDeltaPhiDeg_pfx;
Double_t heading0;
Double_t funcToMin(Double_t pitchOffet, Double_t rollOffset, Double_t headingOffset);


const Int_t numPhiDegBins = 128; //1024;
const Double_t maxPhiDeg = 360;


int main(int argc, char *argv[])
{
  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;    
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  

  TChain* chain = new TChain("angResTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    chain->Add(TString::Format("generateAngularResolutionTree_run%d-%dPlots_noPitchRoll.root", run, run));
  }

  TFile* gpsFile = TFile::Open("~/UCL/ANITA/flight1415/root/run333/gpsFile333.root");
  TTree* adu5PatTree = (TTree*) gpsFile->Get("adu5PatTree");
  Adu5Pat* pat = NULL;
  adu5PatTree->SetBranchAddress("pat", &pat);


  
  TString outFileName = TString::Format("%sPlots.root", argv[0]);
  TFile* outFile = new TFile(outFileName, "recreate");


  TH2D* hDeltaThetaDeg = new TH2D("hDeltaThetaDeg",
  				  "#delta#theta vs #phi_{expected}",
  				  numPhiDegBins, 0, maxPhiDeg,
  				  128, -5, 5);				  
  hDeltaThetaDeg->GetXaxis()->SetTitle("#phi_{expected} (Degrees)");
  hDeltaThetaDeg->GetYaxis()->SetTitle("#delta#theta (Degrees)");
  // chain->Draw("deltaThetaDeg:phiExpected>>hDeltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 5", "goff");
  chain->Draw("deltaThetaDeg:zoomPhiDeg>>hDeltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 5", "goff");  

  hDeltaThetaDeg_pfx = hDeltaThetaDeg->ProfileX();
  hDeltaThetaDeg_pfx->GetYaxis()->SetTitle(hDeltaThetaDeg->GetYaxis()->GetTitle());  

  TH2D* hDeltaPhiDeg = new TH2D("hDeltaPhiDeg",
  				  "#delta#phi vs #phi_{expected}",
  				  numPhiDegBins, 0, maxPhiDeg,
  				  128, -5, 5);				  
  hDeltaPhiDeg->GetXaxis()->SetTitle("#phi_{expected} (Degrees)");
  hDeltaPhiDeg->GetYaxis()->SetTitle("#delta#phi (Degrees)");
  // chain->Draw("deltaPhiDeg:phiExpected>>hDeltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 5", "goff");
  chain->Draw("deltaPhiDeg:zoomPhiDeg>>hDeltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 5", "goff");

  hDeltaPhiDeg_pfx = hDeltaPhiDeg->ProfileX();
  hDeltaPhiDeg_pfx->GetYaxis()->SetTitle(hDeltaPhiDeg->GetYaxis()->GetTitle());  

  Double_t mean=0;
  for(Int_t binx=1; binx<=hDeltaPhiDeg_pfx->GetNbinsX(); binx++){
    mean += hDeltaPhiDeg_pfx->GetBinContent(binx);
    // std::cout << hDeltaPhiDeg_pfx->GetBinContent(binx) << std::endl;
  }
  std::cout << mean/hDeltaPhiDeg_pfx->GetNbinsX() << std::endl;


  adu5PatTree->GetEntry(0);
  usefulPat = UsefulAdu5Pat(pat);
  // usefulPat.pitch=0;
  // usefulPat.roll=0;
  heading0 = usefulPat.heading;
  // Double_t thetaExpected0, phiExpected0;
  // usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0, phiExpected0);

  std::cout << heading0 << "\t" << pat->heading << std::endl;

  // const Double_t minHeadInd = -40;
  // const Double_t maxHeadInd = 40;
  const Double_t minPitchInd = -40;
  const Double_t maxPitchInd = 40;
  const Double_t minRollInd = -40;
  const Double_t maxRollInd = 40;

  Double_t deltaPitch = 0.05;
  Double_t deltaRoll = 0.05;
  // Double_t deltaHead = 0.1;

  Double_t minFom = 1e308;
  Double_t bestPitch = 0;
  // Double_t bestHead = 0;
  Double_t bestRoll = 0;  
  
  Double_t headingOffset = 22; //-0.00664558;
  for(Int_t pitchInd=minPitchInd; pitchInd<maxPitchInd; pitchInd++){
    Double_t pitch = pitchInd * deltaPitch;
    for(Int_t rollInd=minRollInd; rollInd<maxRollInd; rollInd++){
      Double_t roll = rollInd * deltaRoll;

      Double_t fom = funcToMin(pitch, roll, headingOffset);

      if(fom < minFom){
	bestPitch = pitch;
	bestRoll = roll;
	minFom = fom;
      }
    }
  }

  std::cout << "best pitch = " << bestPitch << std::endl;
  std::cout << "best roll = " << bestRoll << std::endl;  

  // bestPitch = 0.05;
  // bestRoll = -1.1;
  
  TGraph* grDeltaTheta = new TGraph();
  for(Int_t headInd=0; headInd < numPhiDegBins; headInd++){

    Double_t thisHeading = headInd*(maxPhiDeg/numPhiDegBins) - heading0;
    if(thisHeading < 0) thisHeading += 360;
    if(thisHeading >= 360) thisHeading -= 360;
    usefulPat.heading = thisHeading;
    usefulPat.pitch = 0;
    usefulPat.roll = 0;

    Double_t thetaExpected0, phiExpected0;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0, phiExpected0);
    thetaExpected0 *= -1*TMath::RadToDeg();
    phiExpected0 *= TMath::RadToDeg();
    

    usefulPat.pitch = bestPitch;
    usefulPat.roll = bestRoll;
    
    Double_t thetaExpected, phiExpected;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
    thetaExpected *= -1*TMath::RadToDeg();
    phiExpected *= TMath::RadToDeg();    
    
    grDeltaTheta->SetPoint(grDeltaTheta->GetN(), phiExpected, thetaExpected0-thetaExpected);
  }
  grDeltaTheta->SetName("grDeltaTheta");
  grDeltaTheta->Write();

  outFile->Write();
  outFile->Close();

  return 0;

}



Double_t funcToMin(Double_t pitchOffset, Double_t rollOffset, Double_t headingOffset){

  Double_t sumOfSquares = 0;
  for(Int_t headInd=0; headInd < numPhiDegBins; headInd++){

    Double_t thisHeading = headInd*(maxPhiDeg/numPhiDegBins) - heading0 + headingOffset;
    usefulPat.heading = thisHeading;

    if(usefulPat.heading < 0) usefulPat.heading += 360;
    if(usefulPat.heading >= 360) usefulPat.heading -= 360;
    usefulPat.pitch = 0;
    usefulPat.roll = 0;

    Double_t thetaExpected0, phiExpected0;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected0, phiExpected0);
    thetaExpected0 *= -1*TMath::RadToDeg();
    phiExpected0 *= TMath::RadToDeg();

    usefulPat.pitch = pitchOffset;
    usefulPat.roll = rollOffset;
    usefulPat.heading = thisHeading + headingOffset;
    if(usefulPat.heading < 0) usefulPat.heading += 360;
    if(usefulPat.heading >= 360) usefulPat.heading -= 360;
    
    Double_t thetaExpected, phiExpected;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
    thetaExpected *= -1*TMath::RadToDeg();
    phiExpected *= TMath::RadToDeg();

    Double_t phiBin = hDeltaThetaDeg_pfx->FindBin(phiExpected);
    Double_t deltaThetaMeasured = hDeltaThetaDeg_pfx->GetBinContent(phiBin); // OldExpected - measured
    Double_t newDeltaTheta = thetaExpected0 - thetaExpected; // OldExpected - newExpected

    // Double_t deltaPhiMeasured = hDeltaPhiDeg_pfx->GetBinContent(phiBin); // OldExpected - measured
    // Double_t newDeltaPhi = phiExpected0 - phiExpected; // OldExpected - newExpected

    // std::cout << deltaThetaMeasured << "\t" << deltaPhiMeasured << std::endl;
    
    // Want newExpected - measured
    // (OldExpected - measured) - (OldExpected - newExpected) = newExpected - measured
    Double_t diffTheta = deltaThetaMeasured - newDeltaTheta;
    // Double_t diffPhi = deltaPhiMeasured - newDeltaPhi;
    Double_t diffSquared = diffTheta*diffTheta;// + diffPhi*diffPhi;

    sumOfSquares += diffSquared;
    // grDeltaTheta->SetPoint(grDeltaTheta->GetN(), phiExpected, thetaExpected0-thetaExpected);
  }

  return sumOfSquares;
}
