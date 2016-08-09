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


double optAllAnt(double *par);
vector<double> findSlope(double **xIn, double **yIn, int *count);
// vector<double> leastSquares(double *xIn, double *yIn, int count);
double getMean(double **x, int ant, int n);
double getRMS(double **x, int ant, int n);

double optAllAnt(double *par){

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
//   Double_t deltaHeading[1]={0};

  for(unsigned int i = 0; i<MAX_ANTENNAS;++i){
    deltaR[i]=par[i];
    deltaZ[i]=par[i+MAX_ANTENNAS];
    deltaPhi[i]=par[i+MAX_ANTENNAS*2];
    deltaCableDelays[i] = par[i+MAX_ANTENNAS*3];
    //    cout << "ant" << i << ": " << deltaR[i] << " " << deltaZ[i] << " " << deltaPhi[i] << " " << deltaCableDelays[i] << endl;

  }

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

  double **AdjecentAntDeltaT;
  double **AdjecentAntPhiWave;
  double **VerticalAntDeltaT;
  double **VerticalAntPhiWave;


  AdjecentAntDeltaT = new double*[48]; // Top ring 0-15, Middle ring 16-31, Bottom ring 32-47
  AdjecentAntPhiWave = new double*[48]; 
  VerticalAntDeltaT = new double*[48]; // Top-Middle ring 0-15, Middle-Bottom ring 16-31, Top-Bottom ring 32-47
  VerticalAntPhiWave = new double*[48]; 
  for (unsigned int i=0;i<48;++i){
    AdjecentAntDeltaT[i] = new double[BIGNUMBER];
    AdjecentAntPhiWave[i] = new double[BIGNUMBER];
    VerticalAntDeltaT[i] = new double[BIGNUMBER];
    VerticalAntPhiWave[i] = new double[BIGNUMBER];
  }

  Int_t countAdjecent[48]={0};
  Int_t countVertical[48]={0};
  Int_t countHorz[48]={0};

  double meanPhi[MAX_ANTENNAS] = {0}; 
  for (unsigned int ant=0;ant<MAX_ANTENNAS;++ant){
    meanPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
    //    meanPhi[ant] = fGeomTool->getAntPhiPosition(ant, pol);
  }

  Double_t additionalPhi = 22.5*TMath::DegToRad();

  Double_t phiWave, thetaWave, deltaTExpected, maxCorrTime;

  Float_t lower, upper;
  Float_t TwoPi = TMath::Pi()*2;

  unsigned int ant1, ant2, vert;

  for(Long64_t entry=0;entry<maxEntry;++entry) {
//     if(entry%100==0) std::cout << entry*100./maxEntry << "%        \r" << flush;
    
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

      //        if (cor->maxCorVals[corInd]/cor->rmsCorVals[corInd]<6) continue;
      maxCorrTime = cor->maxCorTimes[corInd];
      ant1 = cor->firstAnt[corInd];
      ant2 = cor->secondAnt[corInd];

      deltaTExpected=usefulPat.getDeltaTExpectedOpt(ant1, ant2, sourceLon, sourceLat, sourceAlt, deltaR, deltaZ, deltaPhi) + deltaCableDelays[ant1] - deltaCableDelays[ant2];

       
       if (corInd>5 && corInd<12){

	 lower  = meanPhi[ant1] - additionalPhi ;
	 upper = meanPhi[ant2] + additionalPhi ;
	 if (lower<0) lower+=TwoPi;
	 if (upper>TwoPi) upper-=TwoPi;

	 if (lower>upper){
	   if (phiWave<TwoPi*0.5) lower-=TwoPi;
	   else upper+=TwoPi;
	 }
	 
	 if (phiWave>lower && phiWave<upper && (maxCorrTime-deltaTExpected)*(maxCorrTime-deltaTExpected)<1){
	   
	   AdjecentAntDeltaT[ant1][countAdjecent[ant1]] = (maxCorrTime - deltaTExpected);
 	   AdjecentAntPhiWave[ant1][countAdjecent[ant1]] = (phiWave);
	   countAdjecent[ant1]++;
	   
	 }
	 
	 
       }
      


       if (corInd<6 || (corInd>36 && corInd<40) ){
	 
	 if (ant1<16){
	   if (ant2<32) vert = ant1;
	   else vert = ant2;
	 } else {
	   vert = ant1;
	 }

	 lower  = meanPhi[ant1] - additionalPhi ;
	 upper = meanPhi[ant2] + additionalPhi ;
	 if (lower<0) lower+=TwoPi;
	 if (upper>TwoPi) upper-=TwoPi;

	 if (lower>upper){
	   if (phiWave<TwoPi*0.5) lower-=TwoPi;
	   else upper+=TwoPi;
	 }


  	 if ( phiWave>lower && phiWave<upper && (maxCorrTime-deltaTExpected)*(maxCorrTime-deltaTExpected)<1){

	 
	   VerticalAntDeltaT[vert][countVertical[vert]] = (maxCorrTime - deltaTExpected);
 	   VerticalAntPhiWave[vert][countVertical[vert]] = (phiWave);
	   countVertical[vert]++;
	 }
	 
       }

     }

  }

  delete header;
  delete pat;
  delete cor;
  delete gpsChain;
  delete headChain;
  delete corChain;
//   delete fGeomTool;

  Double_t sumMeanAdj = 0;
  Double_t sumMeanVert = 0;
  Double_t sumRMSAdj = 0;
  Double_t sumRMSVert = 0;
  Double_t sumGRADAdj = 0;
  Double_t sumGRADVert = 0;
  double mean = 0;
  double slope = 0;

  vector<double> slopeADJ = findSlope(AdjecentAntPhiWave, AdjecentAntDeltaT, countAdjecent);
  vector<double> slopeVERT = findSlope(VerticalAntPhiWave, VerticalAntDeltaT, countVertical);
  for(unsigned int ants = 0; ants < MAX_ANTENNAS; ++ants){
    
    mean = getMean(AdjecentAntDeltaT, ants, countAdjecent[ants]);
    sumMeanAdj += mean*mean;
//     sumRMSAdj += getRMS(AdjecentAntDeltaT, ants, countAdjecent[ants]);

    sumGRADAdj += slopeADJ[ants]*slopeADJ[ants];

    delete [] AdjecentAntPhiWave[ants];
    delete [] AdjecentAntDeltaT[ants];   

    mean = getMean(VerticalAntDeltaT, ants, countVertical[ants]);
    sumMeanVert += mean*mean;
//     sumRMSVert += getRMS(VerticalAntDeltaT, ants, countVertical[ants]);

//     slope= findSlope(VerticalAntPhiWave, VerticalAntDeltaT, ants, countVertical[ants] );	
    sumGRADVert += slopeVERT[ants]*slopeVERT[ants];

    delete [] VerticalAntPhiWave[ants];
    delete [] VerticalAntDeltaT[ants];    
  }
    
  
  delete [] AdjecentAntDeltaT;
  delete [] AdjecentAntPhiWave;
  delete [] VerticalAntDeltaT;
  delete [] VerticalAntPhiWave;

  
  cout << sumMeanAdj << "  " << sumGRADAdj << " " << sumMeanVert << " " << sumGRADVert << endl;
  
  return (sumMeanAdj+sumGRADAdj+sumMeanVert+sumGRADVert);
   
}

// vector<double> leastSquares(double *xIn, double *yIn, int count) {
  
//   double SUMx, SUMy, SUMxy, SUMxx, SUMres, res, slope,
//     y_intercept, y_estimate ;
//   int i,n;
  
//   SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
//   for (i=0; i<count; i++) {
    
//     SUMx = SUMx + xIn[i];
//     SUMy = SUMy + yIn[i];
//     SUMxy = SUMxy + xIn[i]*yIn[i];
//     SUMxx = SUMxx + xIn[i]*xIn[i];
//   }
  
//   slope = ( SUMx*SUMy - count*SUMxy ) / ( SUMx*SUMx - count*SUMxx );
//   y_intercept = ( SUMy - slope*SUMx ) / count;


//   vector<double> theReturn2;
//   theReturn2.push_back(slope);
//   theReturn2.push_back(y_intercept);
//   return theReturn2;
  
// }


vector<double> findSlope(double **xIn, double **yIn, int *count) {
  
  double SUMx, SUMy, SUMxy, SUMxx, SUMres, slope, y_intercept;
  
  SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
  vector<double> theReturn;

  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    for (unsigned int i=0; i<count[ant]; ++i) {
      
      SUMx = SUMx + xIn[ant][i];
      SUMy = SUMy + yIn[ant][i];
      SUMxy = SUMxy + xIn[ant][i]*yIn[ant][i];
      SUMxx = SUMxx + xIn[ant][i]*xIn[ant][i];
    }
    slope = ( SUMx*SUMy - count[ant]*SUMxy ) / ( SUMx*SUMx - count[ant]*SUMxx );
    theReturn.push_back(slope);
  }


//   cout << slope << endl;

//   y_intercept = ( SUMy - slope*SUMx ) / count;

//   theReturn2.push_back(slope);
//   theReturn2.push_back(y_intercept);
  return theReturn;
  
}




double getMean(double **x, int ant, int n){

  double mean = 0;
  for (int i=0;i<n;++i) mean+=x[ant][i];

  return mean/n;

}


double getRMS(double **x, int ant, int n){

  double rms = 0;
  for (int i=0;i<n;++i) rms+=x[ant][i]*x[ant][i];

  return rms/n;

}
