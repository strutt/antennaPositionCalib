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
TProfile* hDeltaThetaDeg2_pfx;
TProfile* hDeltaPhiDeg2_pfx;

Double_t heading0;
Double_t funcToMin(Double_t pitchOffet, Double_t rollOffset, Double_t headingOffset);


const Int_t numPhiDegBins = 128; //1024;
const Double_t maxPhiDeg = 360;
Int_t firstRun;
Int_t lastRun;
Int_t numCalls;

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
  firstRun = atoi(argv[1]);
  lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  numCalls = 0;

  TChain* chain = new TChain("angResTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    // chain->Add(TString::Format("generateAngularResolutionTree_run%d-%dPlots_noPitchRoll.root", run, run));
    chain->Add(TString::Format("generateAngularResolutionTree_run%d-%dPlots.root", run, run));
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
  chain->Draw("deltaThetaDeg:phiExpected>>hDeltaThetaDeg", "TMath::Abs(deltaThetaDeg) < 5", "goff");
  hDeltaThetaDeg_pfx = hDeltaThetaDeg->ProfileX();
  hDeltaThetaDeg_pfx->GetYaxis()->SetTitle(hDeltaThetaDeg->GetYaxis()->GetTitle());  


  TH2D* hDeltaThetaDeg2 = new TH2D("hDeltaThetaDeg2",
				   "#delta#theta vs #phi_{zoom}",
				   numPhiDegBins, 0, maxPhiDeg,
				   128, -5, 5);				  
  hDeltaThetaDeg2->GetXaxis()->SetTitle("#phi_{zoom} (Degrees)");
  hDeltaThetaDeg2->GetYaxis()->SetTitle("#delta#theta (Degrees)");
  chain->Draw("deltaThetaDeg:zoomPhiDeg>>hDeltaThetaDeg2", "TMath::Abs(deltaThetaDeg) < 5", "goff");  
  hDeltaThetaDeg2_pfx = hDeltaThetaDeg2->ProfileX();
  hDeltaThetaDeg2_pfx->GetYaxis()->SetTitle(hDeltaThetaDeg2->GetYaxis()->GetTitle());  
  
  
  TH2D* hDeltaPhiDeg = new TH2D("hDeltaPhiDeg",
  				  "#delta#phi vs #phi_{expected}",
  				  numPhiDegBins, 0, maxPhiDeg,
  				  128, -5, 5);				  
  hDeltaPhiDeg->GetXaxis()->SetTitle("#phi_{expected} (Degrees)");
  hDeltaPhiDeg->GetYaxis()->SetTitle("#delta#phi (Degrees)");
  chain->Draw("deltaPhiDeg:phiExpected>>hDeltaPhiDeg", "TMath::Abs(deltaPhiDeg) < 5", "goff");
  hDeltaPhiDeg_pfx = hDeltaPhiDeg->ProfileX();
  hDeltaPhiDeg_pfx->GetYaxis()->SetTitle(hDeltaPhiDeg->GetYaxis()->GetTitle());  

  TH2D* hDeltaPhiDeg2 = new TH2D("hDeltaPhiDeg2",
  				  "#delta#phi vs #phi_{zoom}",
  				  numPhiDegBins, 0, maxPhiDeg,
  				  128, -5, 5);				  
  hDeltaPhiDeg2->GetXaxis()->SetTitle("#phi_{zoom} (Degrees)");
  hDeltaPhiDeg2->GetYaxis()->SetTitle("#delta#phi (Degrees)");
  chain->Draw("deltaPhiDeg:zoomPhiDeg>>hDeltaPhiDeg2", "TMath::Abs(deltaPhiDeg) < 5", "goff");
  hDeltaPhiDeg2_pfx = hDeltaPhiDeg2->ProfileX();
  hDeltaPhiDeg2_pfx->GetYaxis()->SetTitle(hDeltaPhiDeg->GetYaxis()->GetTitle());  
  
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

  const Double_t minHeadInd = 20; //-10;
  const Double_t maxHeadInd = 40; //50;
  const Double_t minPitchInd = -100;
  const Double_t maxPitchInd = 0;
  const Double_t minRollInd = -10;
  const Double_t maxRollInd = 10;

  Double_t deltaPitch = 0.01;
  Double_t deltaRoll = 0.01;
  Double_t deltaHeading = 0.01;

  // best pitch = 0.55
  // best roll = -0.9
  

  Double_t minFom = 1e308;
  Double_t bestPitch = 0;
  Double_t bestHeading = 0;
  Double_t bestRoll = 0;  

  for(Int_t headInd=minHeadInd; headInd<maxHeadInd; headInd++){
    Double_t headingOffset = headInd * deltaHeading;  
    for(Int_t pitchInd=minPitchInd; pitchInd<maxPitchInd; pitchInd++){
      Double_t pitch = pitchInd * deltaPitch;
      for(Int_t rollInd=minRollInd; rollInd<maxRollInd; rollInd++){
	Double_t roll = rollInd * deltaRoll;

	Double_t fom = funcToMin(pitch, roll, headingOffset);
	// std::cout << "fom = " << fom << std::endl;
	if(fom < minFom){
	  bestPitch = pitch;
	  bestRoll = roll;
	  bestHeading = headingOffset;
	  minFom = fom;
	}
      }
    }
  }

  std::cout << "best pitch = " << bestPitch << std::endl;
  std::cout << "best roll = " << bestRoll << std::endl;
  std::cout << "best heading = " << bestHeading << std::endl;    

  // bestPitch = 0; //-0.4; //0.05;
  // bestRoll = 0.4;
  // bestHeading = 0;
  
  TGraph* grDeltaTheta = new TGraph();
  TGraph* grDeltaPhi = new TGraph();
  TGraph* grDeltaThetaNew = new TGraph();
  TGraph* grDeltaPhiNew = new TGraph();  
  
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

    Int_t phiBin = hDeltaPhiDeg_pfx->FindBin(phiExpected0);
    Double_t deltaPhiMeasured = hDeltaPhiDeg_pfx->GetBinContent(phiBin);
    Double_t phiMeasured = phiExpected0 + deltaPhiMeasured;
    Double_t deltaThetaMeasured = hDeltaThetaDeg_pfx->GetBinContent(phiBin);
    Double_t thetaMeasured = thetaExpected0 + deltaThetaMeasured;
    
    if(phiMeasured < 0) phiMeasured += 360;
    if(phiMeasured >= 360) phiMeasured -= 360;

    usefulPat.heading = thisHeading + bestHeading;
    if(usefulPat.heading < 0) usefulPat.heading += 360;
    if(usefulPat.heading >= 360) usefulPat.heading -= 360;
    usefulPat.pitch += bestPitch;
    usefulPat.roll += bestRoll;
    
    Double_t thetaExpected, phiExpected;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
    thetaExpected *= -1*TMath::RadToDeg();
    phiExpected *= TMath::RadToDeg();
    
    grDeltaTheta->SetPoint(grDeltaTheta->GetN(), phiMeasured, thetaExpected-thetaExpected0);
    grDeltaPhi->SetPoint(grDeltaPhi->GetN(), phiMeasured, phiExpected-phiExpected0);
    grDeltaThetaNew->SetPoint(grDeltaThetaNew->GetN(), phiMeasured, thetaMeasured-thetaExpected);
    grDeltaPhiNew->SetPoint(grDeltaPhiNew->GetN(), phiMeasured, phiMeasured-phiExpected);
    
  }
  grDeltaTheta->SetName("grDeltaTheta");
  grDeltaTheta->Sort();
  grDeltaTheta->Write();  
  grDeltaPhi->SetName("grDeltaPhi");
  grDeltaPhi->Sort();
  grDeltaPhi->Write();  
  grDeltaThetaNew->SetName("grDeltaThetaNew");
  grDeltaThetaNew->Sort();
  grDeltaThetaNew->Write();  
  grDeltaPhiNew->SetName("grDeltaPhiNew");
  grDeltaPhiNew->Sort();
  grDeltaPhiNew->Write();  

  outFile->Write();
  outFile->Close();

  return 0;

}



// Double_t funcToMin(Double_t pitchOffset, Double_t rollOffset, Double_t headingOffset){

//   Double_t sumOfSquares = 0;
  
//   TChain* headChain = new TChain("headTree");
//   TChain* gpsChain = new TChain("adu5PatTree");
//   const Double_t maxDeltaTriggerTimeNs = 1200;
  
//   TChain* angResChain = new TChain("angResTree");
  
//   for(Int_t run=firstRun; run<=lastRun; run++){
//     TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
//     headChain->Add(fileName);
//     fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
//     gpsChain->Add(fileName);
//     fileName = TString::Format("generateAngularResolutionTree_run%d-%dPlots_noPitchRoll.root", run, run);
//     angResChain->Add(fileName);
//   }

//   RawAnitaHeader* header = NULL;
//   headChain->SetBranchAddress("header", &header);  
//   Adu5Pat* pat = NULL;
//   gpsChain->SetBranchAddress("pat", &pat);
//   gpsChain->BuildIndex("realTime");


//   UInt_t triggerTimeNs;
//   UInt_t triggerTimeNsExpected;
//   Double_t zoomPhiDeg;
//   Double_t zoomThetaDeg;
//   angResChain->SetBranchAddress("triggerTimeNs", &triggerTimeNs);
//   angResChain->SetBranchAddress("triggerTimeNsExpected", &triggerTimeNsExpected);  
//   angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
//   angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);  

//   Long64_t nEntries = headChain->GetEntries();
//   Long64_t maxEntry = 0;
//   Long64_t startEntry = 0;
//   Long64_t entryCount = 0;
//   if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
//   std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
//   std::cout << "This is call number " << numCalls << std::endl;
//   numCalls++;
//   ProgressBar p(maxEntry-startEntry);

//   for(Long64_t entry = startEntry; entry < maxEntry; entry++){
//     headChain->GetEntry(entry);
//     gpsChain->GetEntryWithIndex(header->realTime);
//     if((header->trigType & 1)==1){
//       usefulPat = UsefulAdu5Pat(pat);
//       UInt_t triggerTimeNsExpectedNew = usefulPat.getWaisDivideTriggerTimeNs();
//       UInt_t triggerTimeNsNew = header->triggerTimeNs;
//       Int_t deltaTriggerTimeNs = Int_t(triggerTimeNsNew) - Int_t(triggerTimeNsExpectedNew);
//       if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
	
// 	angResChain->GetEntry(entryCount);
// 	entryCount++;
	
// 	if(triggerTimeNsExpectedNew!=triggerTimeNsExpected && triggerTimeNsNew != triggerTimeNs){
// 	  std::cout << triggerTimeNsExpectedNew<< "\t" << triggerTimeNsExpected << "\t"
// 		    << triggerTimeNsNew << "\t" << triggerTimeNs << std::endl;
// 	  angResChain->Show(entryCount);
// 	  headChain->Show(entry);
// 	  continue;
// 	}

	
// 	Double_t thetaExpectedNew;
// 	Double_t phiExpectedNew;
// 	usefulPat.heading += headingOffset;
// 	usefulPat.pitch += pitchOffset;	
// 	usefulPat.roll += rollOffset;
	
// 	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpectedNew, phiExpectedNew);
// 	phiExpectedNew*=TMath::RadToDeg();
// 	thetaExpectedNew*=-1*TMath::RadToDeg();

// 	Double_t deltaThetaNew = RootTools::getDeltaAngleDeg(thetaExpectedNew, zoomThetaDeg);
// 	Double_t deltaPhiNew = RootTools::getDeltaAngleDeg(phiExpectedNew, zoomPhiDeg);	
// 	sumOfSquares += deltaThetaNew*deltaThetaNew + deltaPhiNew*deltaPhiNew;
	
//       }
//     }
//     p++;
//   }

//   return sumOfSquares/entryCount;
// }




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


    Int_t phiBin = hDeltaPhiDeg_pfx->FindBin(phiExpected0);
    Double_t deltaPhiMeasured = hDeltaPhiDeg_pfx->GetBinContent(phiBin);
    Double_t phiMeasured = phiExpected0 + deltaPhiMeasured;
    if(phiMeasured < 0) phiMeasured += 360;
    if(phiMeasured >= 360) phiMeasured -= 360;    
    phiBin = hDeltaPhiDeg2_pfx->FindBin(phiMeasured);
    
    Double_t deltaThetaMeasured = hDeltaThetaDeg2_pfx->GetBinContent(phiBin); // OldExpected - measured
    Double_t thetaMeasured = thetaExpected0 + deltaThetaMeasured;

    Double_t diffTheta = RootTools::getDeltaAngleDeg(thetaExpected, thetaMeasured);
    Double_t diffPhi = RootTools::getDeltaAngleDeg(phiExpected, phiMeasured);

    Double_t diffSquared = diffTheta*diffTheta + diffPhi*diffPhi;

    sumOfSquares += diffSquared;
    // grDeltaTheta->SetPoint(grDeltaTheta->GetN(), phiExpected, thetaExpected0-thetaExpected);
  }

  return sumOfSquares;
}
