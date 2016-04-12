// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             The macros don't work with ROOT-6. Time to start again...
             This will recycle as far as possible all the code currently used by Linda.
********************************************************************************************************* */

#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAdu5Pat.h"
#include "CorrelationSummaryAnita3.h"
#include "AnitaGeomTool.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

#include "ProgressBar.h"
#include "OutputConvention.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>


double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave,
			    double *deltaR, double *deltaZ, double *deltaPhi);

// double getMean(std::vector<Double_t>* x);
// double getRMS(std::vector<Double_t>* x);
// double findSlope(std::vector<double> *xIn, std::vector<double> *yIn);

double getMean(const std::vector<Double_t>& x);
double getRMS(const std::vector<Double_t>& x);
double findSlope(Int_t wrappedPhi, const std::vector<double>& xIn, const std::vector<double>& yIn);
  
double quickOptAllAntStepsVPOL(const double *par);

void fillArrays();
void makeGraphs(int fittingStep);



// Set of global vectors which we use to store the variables
std::vector<UInt_t> eventNumberIndex;
std::vector<Double_t> thetaWaveIndex;
std::vector<Double_t> phiWaveIndex;
std::vector<Int_t> antIndex1;
std::vector<Int_t> antIndex2;
std::vector<Double_t> maxCorrTimeIndex;
std::vector<Int_t> adjacent;


Double_t lastVertSlope[NUM_SEAVEYS];
Double_t lastVertMean[NUM_SEAVEYS];
Double_t lastVertRms[NUM_SEAVEYS];

Double_t lastAdjSlope[NUM_SEAVEYS];
Double_t lastAdjMean[NUM_SEAVEYS];
Double_t lastAdjRms[NUM_SEAVEYS];

const Double_t meanScaleFactor = 10000;
const Double_t rmsScaleFactor = 100;
const Double_t sumGradScaleFactor = 10000;

int main(int argc, char *argv[])
{
  
  // read the correlation variables from the TTree into the global arrays.
  fillArrays();


  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }  

  // here I initialize the fitter and give it the function I want to minimize...
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.0001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  const int varsPerAnt = 4; // these are r, z, phi, t
  const int numVars = NUM_SEAVEYS*varsPerAnt;

  Double_t zeros[numVars] = {0};
  quickOptAllAntStepsVPOL(zeros);

  
  ROOT::Math::Functor funcToMin(&quickOptAllAntStepsVPOL, numVars);
  Double_t stepSize = 1e-3;
  std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

  // starting point
  std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);

  min->SetFunction(funcToMin);

  // from inside the fitting function, the variables are unwrapped like so...
  // // deltaR[i]=par[i];
  // // deltaZ[i]=par[i+NUM_SEAVEYS];
  // // deltaPhi[i]=par[i+NUM_SEAVEYS*2];
  // // deltaCableDelays[i] = par[i+NUM_SEAVEYS*3];  


  // so now I give the variables the proper names, values and step sizes.
  for(int varTypeInd = 0; varTypeInd < varsPerAnt; varTypeInd++){
    TString varType;
    switch(varTypeInd){
    case 0:
      varType = "r";
      break;
    case 1:
      varType = "phi";
      break;
    case 2:
      varType = "z";
      break;
    case 3:
      varType = "t";
      break;
    }    
    for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
      const Int_t varInd = varTypeInd*NUM_SEAVEYS + ant;
      TString varName = varType + TString::Format("(%d)", ant);
      variables.at(varInd) = 0;
      min->SetVariable(varInd, std::string(varName.Data()), variables[varInd], step[varInd]);
    }
  }

  const int numFittingSteps = 4;
  for(int fittingStep = 0; fittingStep < numFittingSteps; fittingStep++){
    // first, fix all the variables.
    // I will unfix the relevant ones after.
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      min->FixVariable(ant);
      min->FixVariable(NUM_SEAVEYS + ant);
      min->FixVariable(2*NUM_SEAVEYS + ant);
      min->FixVariable(3*NUM_SEAVEYS + ant);	
    }
    
    switch(fittingStep){

    case 0:
      std::cout << "Fitting ts" << std::endl;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	min->ReleaseVariable(3*NUM_SEAVEYS + ant);	
      }
      break;


    case 1:
      std::cout << "Fitting rs, phis" << std::endl;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	min->ReleaseVariable(ant);
	min->ReleaseVariable(NUM_SEAVEYS + ant);
      }      
      break;


    case 2:
      std::cout << "Fitting ts again" << std::endl;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	min->ReleaseVariable(3*NUM_SEAVEYS + ant);	
      }
      break;

      
    case 3:

      std::cout << "Fitting zs" << std::endl;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	min->ReleaseVariable(2*NUM_SEAVEYS + ant);
      }      
      break;
    } 

    // Always fix t0 = 0;
    min->FixVariable(NUM_SEAVEYS*3);

    
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


    // now we store the intermediate numbers in lovely TGraphs
    makeGraphs(fittingStep);

  }

  outFile->Write();
  outFile->Close();  
  
  return 0;
}






void makeGraphs(int fittingStep){
  // making these graphs gets a little busy so I'll put it out the way


  // Antenna number axis
  Double_t antNums[NUM_SEAVEYS];
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    antNums[ant] = ant;
  }


  // means
  TGraph* grAdjMean = new TGraph(NUM_SEAVEYS, antNums, lastAdjMean);
  TString adjMeanName = TString::Format("grAdjMean_%d", fittingStep);
  grAdjMean->SetName(adjMeanName);
  TString adjMeanTitle = TString::Format("Unweighted adjacent pairs after fitting iteration %d; Antenna; Mean (ns)", fittingStep);
  grAdjMean->SetTitle(adjMeanTitle);
  grAdjMean->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grAdjMean->GetY()[ant];
    grAdjMean->GetY()[ant] = y*y*meanScaleFactor;
  }
  adjMeanTitle = TString::Format("Weighted adjacent pairs after fitting iteration %d; Antenna; Mean (ns)", fittingStep);
  grAdjMean->SetTitle(adjMeanTitle);  
  delete grAdjMean;


  
  TGraph* grVertMean = new TGraph(NUM_SEAVEYS, antNums, lastVertMean);
  TString vertMeanName = TString::Format("grVertMean_%d", fittingStep);
  grVertMean->SetName(vertMeanName);
  TString vertMeanTitle = TString::Format("Unweighted vertical pairs after fitting iteration %d; Antenna; Mean (ns)", fittingStep);
  grVertMean->SetTitle(vertMeanTitle);
    
  grVertMean->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grVertMean->GetY()[ant];
    grVertMean->GetY()[ant] = y*y*meanScaleFactor;
  }
  vertMeanTitle = TString::Format("Weighted vertical pairs after fitting iteration %d; Antenna; Mean (ns)", fittingStep);
  grVertMean->SetTitle(vertMeanTitle);  
  
  delete grVertMean;






  // rms
  TGraph* grAdjRms = new TGraph(NUM_SEAVEYS, antNums, lastAdjRms);
  TString adjRmsName = TString::Format("grAdjRms_%d", fittingStep);
  grAdjRms->SetName(adjRmsName);
  TString adjRmsTitle = TString::Format("Unweighted adjacent pairs after fitting iteration %d; Antenna; Rms (ns)", fittingStep);
  grAdjRms->SetTitle(adjRmsTitle);
  grAdjRms->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grAdjRms->GetY()[ant];
    grAdjRms->GetY()[ant] = y*y*rmsScaleFactor;
  }
  adjRmsTitle = TString::Format("Weighted adjacent pairs after fitting iteration %d; Antenna; Rms (ns)", fittingStep);
  grAdjRms->SetTitle(adjRmsTitle);  
  delete grAdjRms;


  
  TGraph* grVertRms = new TGraph(NUM_SEAVEYS, antNums, lastVertRms);
  TString vertRmsName = TString::Format("grVertRms_%d", fittingStep);
  grVertRms->SetName(vertRmsName);
  TString vertRmsTitle = TString::Format("Unweighted vertical pairs after fitting iteration %d; Antenna; Rms (ns)", fittingStep);
  grVertRms->SetTitle(vertRmsTitle);
    
  grVertRms->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grVertRms->GetY()[ant];
    grVertRms->GetY()[ant] = y*y*rmsScaleFactor;
  }
  vertRmsTitle = TString::Format("Weighted vertical pairs after fitting iteration %d; Antenna; Rms (ns)", fittingStep);
  grVertRms->SetTitle(vertRmsTitle);  
  
  delete grVertRms;

  






  // slope
  TGraph* grAdjSlope = new TGraph(NUM_SEAVEYS, antNums, lastAdjSlope);
  TString adjSlopeName = TString::Format("grAdjSlope_%d", fittingStep);
  grAdjSlope->SetName(adjSlopeName);
  TString adjSlopeTitle = TString::Format("Unweighted adjacent pairs after fitting iteration %d; Antenna; Slope (ns)", fittingStep);
  grAdjSlope->SetTitle(adjSlopeTitle);
  grAdjSlope->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grAdjSlope->GetY()[ant];
    grAdjSlope->GetY()[ant] = y*y*sumGradScaleFactor;
  }
  adjSlopeTitle = TString::Format("Weighted adjacent pairs after fitting iteration %d; Antenna; Slope (ns)", fittingStep);
  grAdjSlope->SetTitle(adjSlopeTitle);  
  delete grAdjSlope;


  
  TGraph* grVertSlope = new TGraph(NUM_SEAVEYS, antNums, lastVertSlope);
  TString vertSlopeName = TString::Format("grVertSlope_%d", fittingStep);
  grVertSlope->SetName(vertSlopeName);
  TString vertSlopeTitle = TString::Format("Unweighted vertical pairs after fitting iteration %d; Antenna; Slope (ns)", fittingStep);
  grVertSlope->SetTitle(vertSlopeTitle);
    
  grVertSlope->Write();

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t y = grVertSlope->GetY()[ant];
    grVertSlope->GetY()[ant] = y*y*sumGradScaleFactor;
  }
  vertSlopeTitle = TString::Format("Weighted vertical pairs after fitting iteration %d; Antenna; Slope (ns/rad)", fittingStep);
  grVertSlope->SetTitle(vertSlopeTitle);  
  
  delete grVertSlope;
  
}




void fillArrays(){
  // UInt_t *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex,
  // 		int *antIndex1, int *antIndex2, double *maxCorrTimeIndex,
  // 		bool *adjacent){

  const int maxNumEvents = 770000; //1500000;
  eventNumberIndex.reserve(maxNumEvents);
  thetaWaveIndex.reserve(maxNumEvents);
  phiWaveIndex.reserve(maxNumEvents);
  antIndex1.reserve(maxNumEvents);
  antIndex2.reserve(maxNumEvents);
  maxCorrTimeIndex.reserve(maxNumEvents);
  adjacent.reserve(maxNumEvents);

  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corName[FILENAME_MAX];
  AnitaGeomTool* fGeomTool = AnitaGeomTool::Instance();  
  fGeomTool->useKurtAnita3Numbers(1);

  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  TChain *corChain = new TChain("corTree");

  for (int run=145;run<162;++run){      
    if (run>149 && run<154){
      continue;
    }
    // sprintf(headerName,"/unix/anita3/flight1415/root/run%d/timedHeadFile%d.root",run,run);
    // sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsEvent%d.root",run,run);
    // sprintf(corName, "/unix/anita3/linda/corTrees/corRun_NEW14_VPOL_%d.root", run);
    sprintf(headerName,"~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%d.root",run,run);
    sprintf(gpsName,"~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root",run,run);
    sprintf(corName, "corTrees/corRun_NEW14_VPOL_%d.root", run);

    headChain->Add(headerName);
    gpsChain->Add(gpsName);
    corChain->Add(corName);
    
  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);
  corChain->SetBranchAddress("cor",&cor);

  headChain->BuildIndex("header->eventNumber");
  gpsChain->BuildIndex("eventNumber");

  const int maxEntry=corChain->GetEntries();

  ProgressBar p(maxEntry);
  
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if(pol == AnitaPol::kVertical){
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    timeOffset = + 92.8;
    //    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }
  else{ 
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    // sourceLat = - (79 + (27.94097/60));
    // sourceLon = -(112 + (6.76208/60));
    // sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }




  Double_t antPhi[NUM_SEAVEYS];
  for (int ant=0;ant<NUM_SEAVEYS;++ant){
    antPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
  }
  
  Double_t additionalPhi = 22.5*TMath::DegToRad();

  Int_t countIndex = 0;  
  for(Long64_t entry=0;entry<maxEntry;++entry) {

    corChain->GetEntry(entry);

    Long64_t headEntry = headChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(headEntry < 0 ){
      continue;
    }
    headChain->GetEntry(headEntry);
    
    if(header->run==151 && header->realTime<1418.9397e6){
      continue;
    }

    //    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);
       
    UsefulAdu5Pat usefulPat(pat);
    
    Double_t phiWave, thetaWave;
    usefulPat.getThetaAndPhiWaveLDB(thetaWave, phiWave);

    for(int corInd=0; corInd<NUM_CORRELATIONS_ANITA3; ++corInd) {      
      
      if (corInd>11 && corInd!=37 && corInd!=38 && corInd!=39){
	continue;
      }

      Double_t maxCorrTime = cor->maxCorTimes[corInd];
      Int_t ant1 = cor->firstAnt[corInd];
      Int_t ant2 = cor->secondAnt[corInd];
      
      double zeros[NUM_SEAVEYS];
      for (int i=0;i<NUM_SEAVEYS;i++){
	zeros[i] = 0;
      }
      Double_t deltaTExpected = getDeltaTExpectedOpt(ant1, ant2, thetaWave, phiWave, zeros, zeros, zeros);     

      // double x[1] = {antPhi[ant1]-phiWave};
      // maxCorrTime -= singleDelayFuncMod(x, parSingleDelay);
      // x[0] = antPhi[ant2] - phiWave;
      // maxCorrTime += singleDelayFuncMod(x, parSingleDelay);

      Double_t lower  = antPhi[ant1] - additionalPhi ;
      Double_t upper = antPhi[ant2] + additionalPhi ;
      if (lower<0){
	lower+=TMath::TwoPi();
      }
      if (upper>TMath::TwoPi()){
	upper-=TMath::TwoPi();
      }

      
      if (lower>upper){
	if (phiWave<TMath::TwoPi()*0.5){
	  lower-=TMath::TwoPi();
	}
	else{
	  upper+=TMath::TwoPi();
	}
      }
            
      if(phiWave>lower && phiWave<upper && (maxCorrTime-deltaTExpected)<1) {

	thetaWaveIndex.push_back(thetaWave);
	phiWaveIndex.push_back(phiWave);
	maxCorrTimeIndex.push_back(maxCorrTime);
	eventNumberIndex.push_back(cor->eventNumber);
	antIndex1.push_back(ant1);
	antIndex2.push_back(ant2);
	// thetaWaveIndex[countIndex] = thetaWave;
	// phiWaveIndex[countIndex] = phiWave;
	// maxCorrTimeIndex[countIndex] = maxCorrTime;
	// eventNumberIndex[countIndex]= cor->eventNumber;
	// antIndex1[countIndex] = ant1;
	// antIndex2[countIndex] = ant2;
	if (corInd>5 && corInd<12){
	  // adjacent[countIndex] = true;
	  adjacent.push_back(1);
	}
	else{
	  adjacent.push_back(0);
	}

	countIndex++;
	// if(countIndex >= maxNumEvents){
	//   std::cerr << "You've got an array bounds problem. That's your fault for coding like it's the 1980s."
	// 	    << std::endl;
	//   std::cerr << "The number of events in your sample is " << countIndex
	// 	    << " but your global array length is " << maxNumEvents << std::endl;
	//   return 1;
	// }
      }
    }

    p++;
  }

  delete header;
  delete pat;
  delete cor;
  delete gpsChain;
  delete headChain;
  delete corChain;

  // antIndex1[countIndex] = -999;  
  std::cout << "Vectors filled with " << countIndex << " events." << std::endl;
}
















double quickOptAllAntStepsVPOL(const double *par){
			       // double *thetaWaveIndex,
			       // double *phiWaveIndex,
			       // int *antIndex1, int *antIndex2,
			       // double *maxCorrTimeIndex,
			       // bool *adjacent){


  // Here we use the photogrammetry numbers
  AnitaGeomTool* fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->useKurtAnita3Numbers(1);

  // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

  Double_t deltaR[NUM_SEAVEYS]={0};  
  Double_t deltaZ[NUM_SEAVEYS]={0};  
  Double_t deltaPhi[NUM_SEAVEYS]={0};  
  Double_t deltaCableDelays[NUM_SEAVEYS]={0};  

  for(int i = 0; i<NUM_SEAVEYS;++i){
    deltaR[i]=par[i];
    deltaZ[i]=par[i+NUM_SEAVEYS];
    deltaPhi[i]=par[i+NUM_SEAVEYS*2];
    deltaCableDelays[i] = par[i+NUM_SEAVEYS*3];
  }

  // Double_t adjacentAntDeltaT[NUM_SEAVEYS][BIGNUMBER];
  // Double_t adjacentAntPhiWave[NUM_SEAVEYS][BIGNUMBER]; 
  // Double_t verticalAntDeltaT[NUM_SEAVEYS][BIGNUMBER];
  // Double_t verticalAntPhiWave[NUM_SEAVEYS][BIGNUMBER];
 
  std::vector<Double_t> adjacentAntDeltaT[NUM_SEAVEYS];
  std::vector<Double_t> adjacentAntPhiWave[NUM_SEAVEYS];
  std::vector<Double_t> verticalAntDeltaT[NUM_SEAVEYS];
  std::vector<Double_t> verticalAntPhiWave[NUM_SEAVEYS];
  for(int ant = 0; ant < NUM_SEAVEYS; ++ant){
    // think this is the right size?
    adjacentAntDeltaT[ant].reserve(eventNumberIndex.size());
    adjacentAntPhiWave[ant].reserve(eventNumberIndex.size());
    verticalAntDeltaT[ant].reserve(eventNumberIndex.size());
    verticalAntPhiWave[ant].reserve(eventNumberIndex.size());
  }

  for(UInt_t entry=0; entry < eventNumberIndex.size(); entry++){
    Double_t thetaWave = thetaWaveIndex.at(entry);
    Double_t phiWave = phiWaveIndex.at(entry);
    
    Double_t maxCorrTime = maxCorrTimeIndex.at(entry);
    Int_t ant1 = antIndex1.at(entry);
    Int_t ant2 = antIndex2.at(entry);
    
    Double_t deltaTExpected = getDeltaTExpectedOpt(ant1, ant2, thetaWave, phiWave,
						   deltaR, deltaZ, deltaPhi);
    deltaTExpected += (deltaCableDelays[ant1] - deltaCableDelays[ant2]);
        
    if (adjacent[entry]){      
      // adjacentAntDeltaT[ant1][countAdjacent[ant1]] = (maxCorrTime - deltaTExpected);
      // adjacentAntPhiWave[ant1][countAdjacent[ant1]] = (phiWave);
      // countAdjacent[ant1]++;
      adjacentAntDeltaT[ant1].push_back(maxCorrTime - deltaTExpected);
      adjacentAntPhiWave[ant1].push_back(phiWave);
    }
    else{
      Int_t vert = 0;
      if (ant1<16){
	if (ant2<32){
	  vert = ant1;
	}
	else{
	  vert = ant2;
	}
      }
      else {
	vert = ant1;
      }
      verticalAntDeltaT[vert].push_back(maxCorrTime - deltaTExpected);
      verticalAntPhiWave[vert].push_back(phiWave);
      // countVertical[vert]++;
      // verticalAntDeltaT[vert][countVertical[vert]] = (maxCorrTime - deltaTExpected);
      // verticalAntPhiWave[vert][countVertical[vert]] = (phiWave);
      // countVertical[vert]++;
      
    }    
    entry++;
  }
 
  // Double_t sumMeanAdj = getMean(adjacentAntDeltaT, countAdjacent)*10000;
  // Double_t sumMeanVert = getMean(verticalAntDeltaT, countVertical)*10000;
  // Double_t sumMeanAdj = getMean(adjacentAntDeltaT)*meanScaleFactor;
  // Double_t sumMeanVert = getMean(verticalAntDeltaT)*meanScaleFactor;  

  Double_t sumMeanAdj = 0;
  Double_t sumMeanVert = 0;
  Double_t sumRMSAdj = 0;
  Double_t sumRMSVert = 0;
  Double_t sumGRADAdj = 0;
  Double_t sumGRADVert = 0;  

  Int_t wrappedPhi[NUM_PHI]={1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};    
  for(int ant=0; ant < NUM_SEAVEYS; ant++){


    // means
    lastAdjMean[ant] = getMean(adjacentAntDeltaT[ant]);
    sumMeanAdj += lastAdjMean[ant]*lastAdjMean[ant];

    lastVertMean[ant] = getMean(verticalAntDeltaT[ant]);
    sumMeanVert += lastVertMean[ant]*lastVertMean[ant];


    // rms
    lastAdjRms[ant] = getRMS(adjacentAntDeltaT[ant]);
    sumRMSAdj += lastAdjRms[ant]*lastAdjRms[ant];

    lastVertRms[ant] = getRMS(verticalAntDeltaT[ant]);
    sumRMSVert += lastVertRms[ant]*lastVertRms[ant];


    // slopes
    lastAdjSlope[ant] = findSlope(wrappedPhi[ant], adjacentAntPhiWave[ant], adjacentAntDeltaT[ant]);
    sumGRADAdj += lastAdjSlope[ant]*lastAdjSlope[ant];

    lastVertSlope[ant] = findSlope(wrappedPhi[ant], verticalAntPhiWave[ant], verticalAntDeltaT[ant]);
    sumGRADVert += lastVertSlope[ant]*lastVertSlope[ant];
    
  }

  // scale factors (global variables)
  sumMeanAdj *= meanScaleFactor;
  sumMeanVert *= meanScaleFactor;
  sumRMSAdj *= rmsScaleFactor;
  sumRMSVert *= rmsScaleFactor;  
  sumGRADAdj *= sumGradScaleFactor;
  sumGRADVert *= sumGradScaleFactor;
    
  std::cout << sumMeanAdj/meanScaleFactor << "  " << sumRMSAdj/rmsScaleFactor << " "
  	    << sumGRADAdj/sumGradScaleFactor << " "
  	    << sumMeanVert/meanScaleFactor << " " << sumRMSVert/rmsScaleFactor << " "
  	    << sumGRADVert/sumGradScaleFactor << std::endl;
  
  return (sumMeanAdj+sumRMSAdj+sumGRADAdj+sumMeanVert+sumRMSVert+sumGRADVert);
   
}


// double findSlope(double **xIn, double **yIn, int *count) {
// double findSlope(std::vector<double> *xIn, std::vector<double> *yIn) {  
  
//   Int_t wrappedPhi[16]={1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};  
//   double theReturn =0;

//   for (int ant=0; ant<NUM_SEAVEYS; ++ant){
//     Double_t SUMx = 0;
//     Double_t SUMy = 0;
//     Double_t SUMxy = 0;
//     Double_t SUMxx = 0;

//     const int n = xIn[ant].size();
//     for(int i=0; i < n; ++i){
//       double thisX = xIn[ant].at(i);
      
//       if(wrappedPhi[ant%16]){
//         if(thisX > TMath::Pi()){
//           thisX -= TMath::TwoPi();
// 	}
//       }

//       SUMx = SUMx + thisX; //xIn[ant][i];
//       SUMy = SUMy + yIn[ant].at(i);
//       SUMxy = SUMxy + thisX*yIn[ant].at(i);
//       SUMxx = SUMxx + thisX*thisX;

//     }
//     Double_t slope = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx );

//     theReturn+=slope*slope;
//   }

//   return theReturn;
  
// }


 double findSlope(Int_t wrappedPhi, const std::vector<double>& xIn, const std::vector<double>& yIn){  

   Double_t SUMx = 0;
   Double_t SUMy = 0;
   Double_t SUMxy = 0;
   Double_t SUMxx = 0;
   
   const int n = xIn.size();
   for(int i=0; i < n; ++i){
     double thisX = xIn.at(i);
      
     if(wrappedPhi){
       if(thisX > TMath::Pi()){
	 thisX -= TMath::TwoPi();
       }
     }

     SUMx = SUMx + thisX; //xIn[i];
     SUMy = SUMy + yIn.at(i);
     SUMxy = SUMxy + thisX*yIn.at(i);
     SUMxx = SUMxx + thisX*thisX;

   }
   Double_t slope = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx );

  return slope;  
}


// double getMean(std::vector<Double_t>* x){
  
//   double sumMean =0;
//   for (unsigned int ant=0; ant<NUM_SEAVEYS; ++ant){
//     double mean = 0;
//     const UInt_t n = x[ant].size();
//     for (UInt_t i=0; i<n; ++i){
//       mean+=x[ant].at(i);
//     }
//     // sumMean += mean*mean/(count[ant]*count[ant]);
//     sumMean += mean*mean/(n*n);    
//   }
//   return sumMean;
// }


// double getMean(std::vector<Double_t>* x){
double getMean(const std::vector<Double_t>& x){  

  double mean = 0;
  const UInt_t n = x.size();
  for (UInt_t i=0; i<n; ++i){
    mean+=x.at(i);
  }
  return mean/n;
}


double getRMS(const std::vector<Double_t>& x){

  Double_t rms = 0;
  Double_t mean = 0;    
  const UInt_t n = x.size();
  for (UInt_t i=0; i < n; ++i){
    mean += x.at(i);
    rms += x.at(i)*x.at(i);
  }
  mean /= n;
  rms = rms/n - mean*mean;

  return TMath::Sqrt(rms);

}


// double getRMS(std::vector<Double_t>* x){

//   double sumRms =0;
//   for (int ant=0; ant<NUM_SEAVEYS; ++ant){

//     Double_t rms = 0;
//     Double_t mean = 0;    
//     const UInt_t n = x[ant].size();
//     for (UInt_t i=0; i < n; ++i){
//       mean += x[ant].at(i);
//       rms += x[ant].at(i)*x[ant].at(i);
//     }
//     mean /= n;
//     rms = rms/n - mean*mean;

//     sumRms += rms;
//   }
//   return sumRms;

// }


double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave,
			    double *deltaR, double *deltaZ, double *deltaPhi){
  AnitaGeomTool* fGeomTool = AnitaGeomTool::Instance();  
  //Now fThetaWave and fPhiWave should be correctly set.
  Double_t phi1=fGeomTool->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
  Double_t r1=fGeomTool->getAntR(ant1)+deltaR[ant1];
  Double_t z1=fGeomTool->getAntZ(ant1)+deltaZ[ant1];

  //    Double_t phi1 = antPhi[ant1]+deltaPhi[ant1];
  //    Double_t r1   = antR[ant1]+deltaR[ant1];
  //    Double_t z1   = antZ(ant1)+deltaZ[ant1];

  Double_t phi2=fGeomTool->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
  Double_t r2=fGeomTool->getAntR(ant2)+deltaR[ant2];
  Double_t z2=fGeomTool->getAntZ(ant2)+deltaZ[ant2];

  //   std::std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::std::endl;

  Double_t tanThetaW=TMath::Tan(thetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(phiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(phiWave-phi2);
   
  return  1e9*((TMath::Cos(thetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}
