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

double cableOptAllAnt(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);
double findSlope(double **xIn, double **yIn, int *count);
double getMean(double **x, int *n);
double getRMS(double **x, int *n);

void  fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);

double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi);

double cableOptAllAnt(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent){

  const int BIGNUMBER = 18500;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if (pol == AnitaPol::kVertical){
    sourceLat = - (79 + (27.93728/60));
    sourceLon = -(112 + (6.74974/60));
    sourceAlt = 1813.42;
    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }else{ 
    sourceLat = - (79 + (27.94097/60));
    sourceLon = -(112 + (6.76208/60));
    sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }

  Double_t deltaR[MAX_ANTENNAS]={0};  
  Double_t deltaZ[MAX_ANTENNAS]={0};  
  Double_t deltaPhi[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays[MAX_ANTENNAS]={0};  

  for(unsigned int i = 0; i<MAX_ANTENNAS;++i){
    deltaCableDelays[i] = par[i];
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
 

  Double_t sumMeanAdj = getMean(AdjacentAntDeltaT, countAdjacent)*10000.;
  Double_t sumMeanVert = getMean(VerticalAntDeltaT, countVertical)*10000.;
  Double_t sumRMSAdj = 0;//getRMS(AdjacentAntDeltaT, countAdjacent);
  Double_t sumRMSVert = 0;//getRMS(VerticalAntDeltaT, countVertical);
  Double_t sumGRADAdj = 0;//findSlope(AdjacentAntPhiWave, AdjacentAntDeltaT, countAdjacent);
  Double_t sumGRADVert = 0;//findSlope(VerticalAntPhiWave, VerticalAntDeltaT, countVertical);

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

  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    for (unsigned int i=0; i<count[ant]; ++i) {
      
      SUMx = SUMx + xIn[ant][i];
      SUMy = SUMy + yIn[ant][i];
      SUMxy = SUMxy + xIn[ant][i]*yIn[ant][i];
      SUMxx = SUMxx + xIn[ant][i]*xIn[ant][i];
    }
    slope = ( SUMx*SUMy - count[ant]*SUMxy ) / ( SUMx*SUMx - count[ant]*SUMxx );
//   y_intercept = ( SUMy - slope*SUMx ) / count;
//     theReturn.push_back(slope);
    theReturn+=slope*slope;
  }

  return theReturn*10000.;
  
}
 

double getMean(double **x, int *count){

  double mean = 0;
  double sumMean =0;
  double rms = 0;
  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    mean = rms = 0;
    for (int i=0;i<count[ant];++i){
      mean+=x[ant][i];
      //      rms+=x[ant][i]*x[ant][i];
    }
    //    sumMean += (mean*mean/(count[ant]*count[ant]))/(rms/count[ant]);
    sumMean += (mean*mean/(count[ant]*count[ant]));


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
 char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corName[FILENAME_MAX];

  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  TChain *corChain = new TChain("corTree");

  for (unsigned int run=331;run<356;++run){
      
    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/headFile%d.root",run,run);
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsFile%d.root",run,run);
    sprintf(corName, "/home/lindac/ANITA/Software/EventCorrelator/macros/corTrees/corRun_NEW3_HPOL_%d.root", run);
    
    headChain->Add(headerName);
    gpsChain->Add(gpsName);
    corChain->Add(corName);
    
  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);
  corChain->SetBranchAddress("cor",&cor);

  headChain->BuildIndex("header->eventNumber");
  gpsChain->BuildIndex("pat->realTime");

  int maxEntry=corChain->GetEntries();
  

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if (pol == AnitaPol::kVertical){
    sourceLat = - (79 + (27.93728/60));
    sourceLon = -(112 + (6.74974/60));
    sourceAlt = 1813.42;
    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }else{ 
    sourceLat = - (79 + (27.94097/60));
    sourceLon = -(112 + (6.76208/60));
    sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }


  Double_t phiWave, thetaWave, lower, upper;
  Double_t maxCorrTime, deltaTExpected;
  int ant1, ant2;

  int countIndex = 0;
  int countNotsaved = 0;

  bool save = false;


  for (unsigned int ant=0;ant<MAX_ANTENNAS;++ant){
    antPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
    //    meanPhi[ant] = fGeomTool->getAntPhiPosition(ant, pol);
  }
  Double_t additionalPhi = 22.5*TMath::DegToRad();
  Double_t TwoPi = TMath::Pi()*2.;
  for(Long64_t entry=0;entry<maxEntry;++entry) {
    save=false;
    corChain->GetEntry(entry);

    Long64_t headEntry = headChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(headEntry < 0 ) continue;
    headChain->GetEntry(headEntry);    
    
    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);
       


    UsefulAdu5Pat usefulPat(pat);
    
    usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaWave, phiWave);


    for(unsigned int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;++corInd) {
      
      
      if (corInd>11 && corInd!=37 && corInd!=38 && corInd!=39) continue;

      maxCorrTime = cor->maxCorTimes[corInd];
      ant1 = cor->firstAnt[corInd];
      ant2 = cor->secondAnt[corInd];
      
      deltaTExpected=usefulPat.getDeltaTExpected(ant1, ant2, sourceLon, sourceLat, sourceAlt);

      lower  = antPhi[ant1] - additionalPhi ;
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
	if (corInd>5 && corInd<12) adjacent[countIndex] = true;
	else adjacent[countIndex]=false;

	countIndex++;
	
      }
    }

	
  }

  delete header;
  delete pat;
  delete cor;
  delete gpsChain;
  delete headChain;
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

//    Double_t phi1 = antPhi[ant1]+deltaPhi[ant1];
//    Double_t r1   = antR[ant1]+deltaR[ant1];
//    Double_t z1   = antZ(ant1)+deltaZ[ant1];

   Double_t phi2=fGeomTool->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
   Double_t r2=fGeomTool->getAntR(ant2)+deltaR[ant2];
   Double_t z2=fGeomTool->getAntZ(ant2)+deltaZ[ant2];

   //   std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::endl;

   Double_t tanThetaW=TMath::Tan(thetaWave);
   Double_t part1=z1*tanThetaW - r1 * TMath::Cos(phiWave-phi1);
   Double_t part2=z2*tanThetaW - r2 * TMath::Cos(phiWave-phi2);
   
   return  1e9*((TMath::Cos(thetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}
