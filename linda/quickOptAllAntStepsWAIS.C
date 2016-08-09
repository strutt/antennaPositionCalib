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

#define MAX_ANTENNAS 48

AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
double antPhi[MAX_ANTENNAS] = {0}; 

double quickOptAllAntStepsWAIS(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);
double findSlope(double **xIn, double **yIn, int *count);
double getMean(double **x, int *n);
double getRMS(double **x, int *n);
Double_t singleDelayFuncMod(Double_t *x, Double_t *par);

void  fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);

double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi);

double quickOptAllAntStepsWAIS(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent){

  fGeomTool->useKurtAnita3Numbers(1);

  const int BIGNUMBER = 100000;//18500;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  // Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if (pol == AnitaPol::kVertical){
    // sourceLat = - (79 + (27.93728/60));
    // sourceLon = -(112 + (6.74974/60));
    // sourceAlt = 1813.42;
    // timeOffset = + 92.8;
    //    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }else{ 
    // sourceLat = - (79 + (27.94097/60));
    // sourceLon = -(112 + (6.76208/60));
    // sourceAlt = 1819.62;
    // timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }

  Double_t deltaR[MAX_ANTENNAS]={0};  
  Double_t deltaZ[MAX_ANTENNAS]={0};  
  Double_t deltaPhi[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays[MAX_ANTENNAS]={0};  

  for(unsigned int i = 0; i<MAX_ANTENNAS;++i){
    deltaR[i]=par[i];
    deltaZ[i]=par[i+MAX_ANTENNAS];
    deltaPhi[i]=par[i+MAX_ANTENNAS*2];
    deltaCableDelays[i] = par[i+MAX_ANTENNAS*3];
    //    cout << "ant" << i << ": " << deltaR[i] << " " << deltaZ[i] << " " << deltaPhi[i] << " " << deltaCableDelays[i] << endl;

  }

  double **AdjacentAntDeltaT;
  double **AdjacentAntPhiWave;
  double **VerticalAntDeltaT;
  double **VerticalAntPhiWave;


  AdjacentAntDeltaT = new double*[48]; // Top ring 0-15, Middle ring 16-31, Bottom ring 32-47
  AdjacentAntPhiWave = new double*[48]; 
  VerticalAntDeltaT = new double*[48]; // Top-Middle ring 0-15, Middle-Bottom ring 16-31, Top-Bottom ring 32-47
  VerticalAntPhiWave = new double*[48]; 
  for (unsigned int i=0;i<48;++i){
    AdjacentAntDeltaT[i] = new double[BIGNUMBER];
    AdjacentAntPhiWave[i] = new double[BIGNUMBER];
    VerticalAntDeltaT[i] = new double[BIGNUMBER];
    VerticalAntPhiWave[i] = new double[BIGNUMBER];
  }

  Int_t countAdjacent[48]={0};
  Int_t countVertical[48]={0};
  Int_t countHorz[48]={0};

  for (unsigned int ant=0;ant<MAX_ANTENNAS;++ant){
    antPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
    //    meanPhi[ant] = fGeomTool->getAntPhiPosition(ant, pol);
  }

  Double_t additionalPhi = 22.5*TMath::DegToRad();

  Double_t phiWave, thetaWave, deltaTExpected, maxCorrTime;

  Float_t lower, upper;
  Float_t TwoPi = TMath::Pi()*2;

  unsigned int ant1, ant2, vert;

  int arrayCount =0;


  Int_t entry =0;

  while(antIndex1[entry]!=-999) {
    //    std::cout << entry << "%        \r" << flush;    
    
    thetaWave = thetaWaveIndex[entry];
    phiWave = phiWaveIndex[entry];
    
    maxCorrTime = maxCorrTimeIndex[entry];
    ant1 = antIndex1[entry];
    ant2 = antIndex2[entry];
    
    //       deltaTExpected=usefulPat.getDeltaTExpectedOpt(ant1, ant2, sourceLon, sourceLat, sourceAlt, deltaR, deltaZ, deltaPhi) + deltaCableDelays[ant1] - deltaCableDelays[ant2];
    
    deltaTExpected = getDeltaTExpectedOpt(ant1, ant2, thetaWave, phiWave, deltaR, deltaZ, deltaPhi) + deltaCableDelays[ant1] - deltaCableDelays[ant2];
    
//     if ((maxCorrTime-deltaTExpected)*(maxCorrTime-deltaTExpected)>1){
//       entry++;
//       continue;
//     }
    //     lower  = antPhi[ant1] - additionalPhi ;
    //     upper = antPhi[ant2] + additionalPhi ;
    //     if (lower<0) lower+=TwoPi;
    //     if (upper>TwoPi) upper-=TwoPi;
    
    //     if (lower>upper){
    //       if (phiWave<TwoPi*0.5) lower-=TwoPi;
    //       else upper+=TwoPi;
    //     }
    
    //     if (phiWave<lower || phiWave>upper){
    //       entry++;
    //       continue;
    //     }
    
    if (adjacent[entry]){
      
      AdjacentAntDeltaT[ant1][countAdjacent[ant1]] = (maxCorrTime - deltaTExpected);
      AdjacentAntPhiWave[ant1][countAdjacent[ant1]] = (phiWave);
      countAdjacent[ant1]++;
    }else {
      
      if (ant1<16){
	if (ant2<32) vert = ant1;
	else vert = ant2;
      } else {
	vert = ant1;
      }
      VerticalAntDeltaT[vert][countVertical[vert]] = (maxCorrTime - deltaTExpected);
      VerticalAntPhiWave[vert][countVertical[vert]] = (phiWave);
      countVertical[vert]++; 
    }
    
    entry++;
  }
 

  Double_t sumMeanAdj = getMean(AdjacentAntDeltaT, countAdjacent)*10000;
  Double_t sumMeanVert = getMean(VerticalAntDeltaT, countVertical)*10000;
  Double_t sumRMSAdj = getRMS(AdjacentAntDeltaT, countAdjacent)*100;
  Double_t sumRMSVert = getRMS(VerticalAntDeltaT, countVertical)*100;
  Double_t sumGRADAdj = findSlope(AdjacentAntPhiWave, AdjacentAntDeltaT, countAdjacent)*10000;
  Double_t sumGRADVert = findSlope(VerticalAntPhiWave, VerticalAntDeltaT, countVertical)*10000;

  for(unsigned int ants = 0; ants < MAX_ANTENNAS; ++ants){

    delete [] AdjacentAntPhiWave[ants];
    delete [] AdjacentAntDeltaT[ants];   

    delete [] VerticalAntPhiWave[ants];
    delete [] VerticalAntDeltaT[ants];    
  }
    
  
  delete [] AdjacentAntDeltaT;
  delete [] AdjacentAntPhiWave;
  delete [] VerticalAntDeltaT;
  delete [] VerticalAntPhiWave;

  
  cout << sumMeanAdj << "  " << sumRMSAdj << " " << sumGRADAdj << " " << sumMeanVert << " " << sumRMSVert << " " << sumGRADVert << endl;
  
  return (sumMeanAdj+sumRMSAdj+sumGRADAdj+sumMeanVert+sumRMSVert+sumGRADVert);
   
}

double findSlope(double **xIn, double **yIn, int *count) {
  
  double SUMx, SUMy, SUMxy, SUMxx, SUMres, slope, y_intercept;
  
  SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
  double theReturn =0;

  Int_t wrappedPhi[16]={1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
    for (unsigned int i=0; i<count[ant]; ++i) {
      double thisX=xIn[ant][i];

      if (wrappedPhi[ant%16]){
	if (thisX>TMath::Pi())
	  thisX-=TMath::TwoPi();
      }

      SUMx = SUMx + thisX; //xIn[ant][i];
      SUMy = SUMy + yIn[ant][i];
      SUMxy = SUMxy + thisX*yIn[ant][i];
      SUMxx = SUMxx + thisX*thisX;
    }
    slope = ( SUMx*SUMy - count[ant]*SUMxy ) / ( SUMx*SUMx - count[ant]*SUMxx );
//   y_intercept = ( SUMy - slope*SUMx ) / count;
//     theReturn.push_back(slope);
    theReturn+=slope*slope;
  }

  return theReturn;
  
}
 



double getMean(double **x, int *count){

  double mean = 0;
  double sumMean =0;
  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    mean = 0;
    for (int i=0;i<count[ant];++i) mean+=x[ant][i];
    sumMean += mean*mean/(count[ant]*count[ant]);
  }
  return sumMean;

}


double getRMS(double **x, int *count){

  double rms = 0;
  double sumRms =0;
  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    rms = 0;
    for (int i=0;i<count[ant];++i) rms+=x[ant][i]*x[ant][i];
    sumRms += rms/count[ant];
  }
  return sumRms;

}


void fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent){

  char gpsName[FILENAME_MAX];
  char corName[FILENAME_MAX];
  // char pointName[FILENAME_MAX];
  fGeomTool->useKurtAnita3Numbers(1);

  Adu5Pat *pat =0;
  CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *corChain = new TChain("corTree");
  // TChain *pointChain = new TChain("pointed");


  for (unsigned int run=331;run<356;++run){
      
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsEvent%d.root",run,run);
    sprintf(corName, "/unix/anita3/linda/corTrees/corRun_NEW12_HPOL_%d.root", run);
    // sprintf(pointName, "/unix/anita3/linda/angularResolution/TreesFromCosmin/out_WAIS_022916/%d.root", run);
    gpsChain->Add(gpsName);
    corChain->Add(corName);
    // pointChain->Add(pointName);
  }
  gpsChain->SetBranchAddress("pat",&pat);
  corChain->SetBranchAddress("cor",&cor);

  // pointChain->SetMakeClass(1);

  // Double_t snrCoherent;
  // Int_t peventNumber;
  // pointChain->SetBranchAddress("eventNumber", &peventNumber);
  // pointChain->SetBranchAddress("snrCoherent", &snrCoherent);

  // pointChain->BuildIndex("eventNumber");
  gpsChain->BuildIndex("eventNumber");

  int maxEntry=corChain->GetEntries();
  

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if (pol == AnitaPol::kVertical){
    sprintf(cpol, "VPOL");
  }else{ 
    sprintf(cpol, "HPOL");
  }


  Double_t phiWave, thetaWave, lower, upper;
  Double_t maxCorrTime, deltaTExpected;
  int ant1, ant2;

  int countIndex = 0;
  int countNotsaved = 0;

  bool save = false;


  Double_t parSingleDelay[13] ={0.00000e+00, 
				0.00000e+00 ,
			       -1.58903e-05 ,
				0.00000e+00 ,
				2.57782e-08 ,
				0.00000e+00 ,
			       -6.98967e-12 ,
				0.00000e+00 ,
				7.98384e-16 ,
				0.00000e+00 ,
			       -3.60950e-20 ,
				0.00000e+00 ,
				2.94428e-25 }; 

  for (unsigned int ant=0;ant<MAX_ANTENNAS;++ant){
    antPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
  }
  Double_t additionalPhi = 22.5*TMath::DegToRad();//22.5*TMath::DegToRad();
  Double_t TwoPi = TMath::Pi()*2.;
  for(Long64_t entry=0;entry<maxEntry;++entry) {
    save=false;
    corChain->GetEntry(entry);

    // Long64_t pointEntry = pointChain->GetEntryNumberWithIndex(cor->eventNumber);
    // if(pointEntry < 0 ) continue;
    // pointChain->GetEntry(pointEntry);
    // if (snrCoherent<5) continue;
    // cout << snrCoherent<< endl;

    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);

    UsefulAdu5Pat usefulPat(pat);
    
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);

    for(unsigned int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;++corInd) {
      
      
      if (corInd>11 && corInd!=37 && corInd!=38 && corInd!=39) continue;
      //      if ( (corInd>11 && corInd<20) || (corInd>25 && corInd<36) && corInd>40) continue;

      // if (cor->maxCorVals[corInd]/cor->rmsCorVals[corInd]<8) continue;

      maxCorrTime = cor->maxCorTimes[corInd];
      ant1 = cor->firstAnt[corInd];
      ant2 = cor->secondAnt[corInd];
      
      deltaTExpected=usefulPat.getDeltaTExpected(ant1, ant2, AnitaLocations::LONGITUDE_WAIS,AnitaLocations::LATITUDE_WAIS,AnitaLocations::ALTITUDE_WAIS);


      double x[1] = {antPhi[ant1]-phiWave};
      maxCorrTime -= singleDelayFuncMod(x, parSingleDelay);
      x[0] = antPhi[ant2] - phiWave;
      maxCorrTime += singleDelayFuncMod(x, parSingleDelay);

      lower = antPhi[ant1] - additionalPhi ;
      upper = antPhi[ant2] + additionalPhi ;
      if (lower<0) lower+=TwoPi;
      if (upper>TwoPi) upper-=TwoPi;

      if (lower>upper){
      	if (phiWave<TwoPi*0.5) lower-=TwoPi;
      	else upper+=TwoPi;
      }
      
      
      if ( phiWave>lower && phiWave<upper && (maxCorrTime-deltaTExpected)*(maxCorrTime-deltaTExpected)<1) {

	thetaWaveIndex[countIndex] = thetaWave;
	phiWaveIndex[countIndex] = phiWave;
	maxCorrTimeIndex[countIndex] = maxCorrTime;
	eventNumberIndex[countIndex]= cor->eventNumber;
	antIndex1[countIndex] = ant1;
	antIndex2[countIndex] = ant2;
	if ( (corInd>5 && corInd<12) || (corInd>19 && corInd<31) ) adjacent[countIndex] = true;
	else adjacent[countIndex]=false;

	countIndex++;
      }
    }

	
  }

  delete pat;
  delete cor;
  delete gpsChain;
  delete corChain;


  antIndex1[countIndex] = -999;

  cout << "ARRAY FILLED! : " << countIndex << " " << countNotsaved << endl;


}



double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi)
{
   //Now fThetaWave and fPhiWave should be correctly set.
    Double_t phi1=fGeomTool->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
    Double_t r1=fGeomTool->getAntR(ant1)+deltaR[ant1];
    Double_t z1=fGeomTool->getAntZ(ant1)+deltaZ[ant1];

   Double_t phi2=fGeomTool->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
   Double_t r2=fGeomTool->getAntR(ant2)+deltaR[ant2];
   Double_t z2=fGeomTool->getAntZ(ant2)+deltaZ[ant2];

   //   std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::endl;

   Double_t tanThetaW=TMath::Tan(thetaWave);
   Double_t part1=z1*tanThetaW - r1 * TMath::Cos(phiWave-phi1);
   Double_t part2=z2*tanThetaW - r2 * TMath::Cos(phiWave-phi2);
   
   return  1e9*((TMath::Cos(thetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}

Double_t singleDelayFuncMod(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=par[1]*TMath::Power(antPhi, 1) +
    + par[2]*TMath::Power(antPhi, 2) 
    + par[3]*TMath::Power(antPhi, 3) 
    + par[4]*TMath::Power(antPhi, 4) 
    + par[5]*TMath::Power(antPhi, 5) 
    + par[6]*TMath::Power(antPhi, 6) 
    + par[7]*TMath::Power(antPhi, 7) 
    + par[8]*TMath::Power(antPhi, 8) 
    + par[9]*TMath::Power(antPhi, 9) 
    + par[10]*TMath::Power(antPhi, 10) 
    + par[11]*TMath::Power(antPhi, 11) 
    + par[12]*TMath::Power(antPhi, 12) 
    // + par[13]*TMath::Power(antPhi, 13) 
    // + par[14]*TMath::Power(antPhi, 14) 
    // + par[15]*TMath::Power(antPhi, 15) 
    // + par[16]*TMath::Power(antPhi, 16) 
    // + par[17]*TMath::Power(antPhi, 17) 
    // + par[18]*TMath::Power(antPhi, 18) 
    // + par[19]*TMath::Power(antPhi, 19) 
    // + par[20]*TMath::Power(antPhi, 20) 
    ;

  return antDelay;
}
