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

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>


#define MAX_ANTENNAS 48

double eventNumberIndex[1500000];
double thetaWaveIndex[1500000];
double phiWaveIndex[1500000];
int antIndex1[1500000];
int antIndex2[1500000];
double maxCorrTimeIndex[1500000];
bool adjacent[1500000];

AnitaGeomTool *fGeomTool = NULL;
double antPhi[MAX_ANTENNAS] = {0};

double quickOptAllAntStepsVPOL(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);
double findSlope(double **xIn, double **yIn, int *count);
double getMean(double **x, int *n);
double getRMS(double **x, int *n);

void  fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent);
Double_t singleDelayFuncMod(Double_t *x, Double_t *par) ;
double getDeltaTExpectedOpt(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi);


void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = quickOptAllAntStepsVPOL(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}



double quickOptAllAntStepsVPOL(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent){
  fGeomTool->useKurtAnita3Numbers(1);
  const int BIGNUMBER = 20000;

  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

  Double_t deltaR[MAX_ANTENNAS]={0};  
  Double_t deltaZ[MAX_ANTENNAS]={0};  
  Double_t deltaPhi[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays[MAX_ANTENNAS]={0};  

  for(unsigned int i = 0; i<MAX_ANTENNAS;++i){
    deltaR[i]=par[i];
    deltaZ[i]=par[i+MAX_ANTENNAS];
    deltaPhi[i]=par[i+MAX_ANTENNAS*2];
    deltaCableDelays[i] = par[i+MAX_ANTENNAS*3];
    //    std::cout << "ant" << i << ": " << deltaR[i] << " " << deltaZ[i] << " " << deltaPhi[i] << " " << deltaCableDelays[i] << std::endl;

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
  // Int_t countHorz[48]={0};

  for (unsigned int ant=0;ant<MAX_ANTENNAS;++ant){
    antPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
    //    meanPhi[ant] = fGeomTool->getAntPhiPosition(ant, pol);
  }

  // Double_t additionalPhi = 22.5*TMath::DegToRad();

  Double_t phiWave, thetaWave, deltaTExpected, maxCorrTime;

  // Float_t lower, upper;
  // Float_t TwoPi = TMath::Pi()*2;

  unsigned int ant1, ant2, vert;

  // int arrayCount =0;


  Int_t entry =0;

  while(antIndex1[entry]!=-999) {
    //    std::std::cout << entry << "         " <<   antIndex1[entry]<< "      \r" << flush;    
    
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
  Double_t sumRMSAdj = 0;//getRMS(AdjacentAntDeltaT, countAdjacent)*100;
  Double_t sumRMSVert = 0;//getRMS(VerticalAntDeltaT, countVertical)*100;
  Double_t sumGRADAdj = 0;//findSlope(AdjacentAntPhiWave, AdjacentAntDeltaT, countAdjacent)*10000;
  Double_t sumGRADVert = 0;//findSlope(VerticalAntPhiWave, VerticalAntDeltaT, countVertical)*10000;

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

  
  std::cout << sumMeanAdj << "  " << sumRMSAdj << " " << sumGRADAdj << " " << sumMeanVert << " " << sumRMSVert << " " << sumGRADVert << std::endl;
  
  return (sumMeanAdj+sumRMSAdj+sumGRADAdj+sumMeanVert+sumRMSVert+sumGRADVert);
   
}

double findSlope(double **xIn, double **yIn, int *count) {
  
  // double SUMx, SUMy, SUMxy, SUMxx, SUMres, slope, y_intercept;
  double SUMx, SUMy, SUMxy, SUMxx, slope;  
  Int_t wrappedPhi[16]={1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};  
  double theReturn =0;

  for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
    SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;

    for (int i=0; i<count[ant]; ++i) {
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
 char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corName[FILENAME_MAX];
  fGeomTool->useKurtAnita3Numbers(1);

  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  TChain *corChain = new TChain("corTree");

  for (unsigned int run=145;run<162;++run){
      
    if (run>149 && run<154) continue;
    //    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/headFile%d.root",run,run);
    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/timedHeadFile%d.root",run,run);
    //    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsFile%d.root",run,run);
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsEvent%d.root",run,run);
    sprintf(corName, "/unix/anita3/linda/corTrees/corRun_NEW14_VPOL_%d.root", run);

    headChain->Add(headerName);
    gpsChain->Add(gpsName);
    corChain->Add(corName);
    
  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);
  corChain->SetBranchAddress("cor",&cor);

  headChain->BuildIndex("header->eventNumber");
  //  gpsChain->BuildIndex("pat->realTime");
  gpsChain->BuildIndex("eventNumber");

  int maxEntry=corChain->GetEntries();
  

  // AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  if (pol == AnitaPol::kVertical){
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    timeOffset = + 92.8;
    //    timeOffset = -99756.6;
    sprintf(cpol, "VPOL");
  }else{ 
    sourceLat = - (77 + (51.23017/60));
    sourceLon = +(167 + (12.16908/60));
    sourceAlt = 0;
    // sourceLat = - (79 + (27.94097/60));
    // sourceLon = -(112 + (6.76208/60));
    // sourceAlt = 1819.62;
    timeOffset = + 92.8;
    sprintf(cpol, "HPOL");
  }


  Double_t phiWave, thetaWave, lower, upper;
  Double_t maxCorrTime, deltaTExpected;
  int ant1, ant2;

  int countIndex = 0;
  int countNotsaved = 0;

  bool save = false;

  double temp[MAX_ANTENNAS];
  for (int i=0;i<MAX_ANTENNAS;i++) temp[i]= 0.;


  // Double_t parSingleDelay[13] ={0.00000e+00, 
  // 				0.00000e+00 ,
  // 				-1.68751e-05 ,
  // 				0.00000e+00 ,
  // 				2.77815e-08 ,
  // 				0.00000e+00 ,
  // 				-8.29351e-12 ,
  // 				0.00000e+00 ,
  // 				1.15064e-15 ,
  // 				0.00000e+00 ,
  // 				-7.71170e-20 ,
  // 				0.00000e+00 ,
  // 				1.99661e-24 }; 


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
    
    if (header->run==151 && header->realTime<1418.9397e6) continue;

    //    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);
       
   // TF1 *fitDelay = new TF1("fitDelay",singleDelayFuncMod,-50,50,5);
   // fitDelay->SetParLimits(0,-1,-1);
   // fitDelay->SetParLimits(1,-1,-1);
   // fitDelay->SetParLimits(3,-1,-1);
   // fitDelay->SetParLimits(5,-1,-1);
   // fitDelay->SetParLimits(7,-1,-1);
   // fitDelay->SetParLimits(9,-1,-1);
   // fitDelay->SetParLimits(11,-1,-1);


    UsefulAdu5Pat usefulPat(pat);
    
    usefulPat.getThetaAndPhiWaveLDB(thetaWave, phiWave);


    for(unsigned int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;++corInd) {
      
      
      if (corInd>11 && corInd!=37 && corInd!=38 && corInd!=39) continue;

      maxCorrTime = cor->maxCorTimes[corInd];
      ant1 = cor->firstAnt[corInd];
      ant2 = cor->secondAnt[corInd];
      
      //      deltaTExpected=usefulPat.getDeltaTExpected(ant1, ant2, sourceLon, sourceLat, sourceAlt);

      deltaTExpected = getDeltaTExpectedOpt(ant1, ant2, thetaWave, phiWave, temp, temp, temp);

      // double x[1] = {antPhi[ant1]-phiWave};
      // maxCorrTime -= singleDelayFuncMod(x, parSingleDelay);
      // x[0] = antPhi[ant2] - phiWave;
      // maxCorrTime += singleDelayFuncMod(x, parSingleDelay);

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

  std::cout << "ARRAY FILLED! : " << countIndex << " " << countNotsaved << std::endl;


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

   //   std::std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::std::endl;

   Double_t tanThetaW=TMath::Tan(thetaWave);
   Double_t part1=z1*tanThetaW - r1 * TMath::Cos(phiWave-phi1);
   Double_t part2=z2*tanThetaW - r2 * TMath::Cos(phiWave-phi2);
   
   return  1e9*((TMath::Cos(thetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


Double_t singleDelayFuncMod(Double_t *x, Double_t *par) {
  Double_t theAntPhi=x[0];
  Double_t antDelay=par[1]*TMath::Power(theAntPhi, 1) +
    + par[2]*TMath::Power(theAntPhi, 2) 
    + par[3]*TMath::Power(theAntPhi, 3) 
    + par[4]*TMath::Power(theAntPhi, 4) 
    + par[5]*TMath::Power(theAntPhi, 5) 
    + par[6]*TMath::Power(theAntPhi, 6) 
    + par[7]*TMath::Power(theAntPhi, 7) 
    + par[8]*TMath::Power(theAntPhi, 8) 
    + par[9]*TMath::Power(theAntPhi, 9) 
    + par[10]*TMath::Power(theAntPhi, 10) 
    + par[11]*TMath::Power(theAntPhi, 11) 
    + par[12]*TMath::Power(theAntPhi, 12) 
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



//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48

int main(){
  
  fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->useKurtAnita3Numbers(1);

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // // set tolerance , etc...
  // min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  // min->SetMaxIterations(10000);  // for GSL
  // min->SetTolerance(0.0001);
  // min->SetPrintLevel(1);

  // // create funciton wrapper for minmizer
  // // a IMultiGenFunction type 
  // // ROOT::Math::Functor funcToMin(&sumOverSquaredDifferences, numVars);
  // ROOT::Math::Functor funcToMin(&quickOptAllAntStepsVPOL, numVars);  

  // Double_t stepSize = 1e-3;
  // std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

  // // starting point
  // std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);

  // min->SetFunction(funcToMin);

  
  TMinuit *myMin = new TMinuit(192);
  myMin->SetObjectFit(quickOptAllAntStepsVPOL);
  myMin->SetFCN(iHateRoot);
 //   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

  double minValCableDelays[MAX_ANTENNAS] ={0};
  double maxValCableDelays[MAX_ANTENNAS] ={0};

  Double_t deltaR[MAX_ANTENNAS] = {0};
  Double_t deltaRErr[MAX_ANTENNAS] = {0};
  Double_t deltaZ[MAX_ANTENNAS] = {0};
  Double_t deltaZErr[MAX_ANTENNAS] = {0};
  Double_t deltaPhi[MAX_ANTENNAS] = {0};
  Double_t deltaPhiErr[MAX_ANTENNAS] = {0};
  Double_t deltaCableDelays[MAX_ANTENNAS] = {0};
  Double_t deltaCableDelaysErr[MAX_ANTENNAS] = {0};

  Double_t globalDeltaz = 0;//-0.0110288; // flip sign for top antenna, same sign for middle bottom antennas
  Double_t globalDeltat = 0.2; // (only for middle and bottom rings)
  
  for(int y = 0; y <MAX_ANTENNAS; y++){
 
    deltaR[y] = 0;
    // minValCableDelays[y] = -0.15 + (y>15)*globalDeltat;
    // maxValCableDelays[y] = +0.15 + (y>15)*globalDeltat;
    deltaZ[y] = (y<16)*(-globalDeltaz) + (y>15)*globalDeltaz;

    // if (y==47){
      minValCableDelays[y] = -0.5;
      maxValCableDelays[y] = +0.5;
    // }
    deltaCableDelays[y] = (y>15)*globalDeltat;
    

    char name[30];
    sprintf(name,"r%d",y);
    myMin->DefineParameter(y, name, deltaR[y], stepSize, -0.3, 0.3);
    sprintf(name,"z%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
    sprintf(name,"phi%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
    sprintf(name,"cable%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
    myMin->FixParameter(y); // fixed R 
    myMin->FixParameter(y+MAX_ANTENNAS); // fixed Z
    myMin->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
    //    myMin->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  myMin->FixParameter(0+MAX_ANTENNAS*3); // fixed t0

  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  


  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    std::cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << std::endl;

    //    myMin->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    //    std::cout << " deltaZ[" << u << "] = " << deltaZ[u] << " +/- " << deltaZErr[u] << std::endl;

  }

  std::cout << " ################## first step done ##############" << std::endl;

  TMinuit *myMin2 = new TMinuit(192);
  myMin2->SetObjectFit(quickOptAllAntStepsVPOL);
  myMin2->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin2->DefineParameter(y, name, deltaR[y], stepSize, -0.3, 0.3);
    sprintf(name,"z%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
    sprintf(name,"phi%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
    sprintf(name,"cable%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
    //    myMin2->FixParameter(y); // fixed R 
    myMin2->FixParameter(y+MAX_ANTENNAS); // fixed Z
    // myMin2->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
    myMin2->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  myMin2->Migrad();


  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin2->GetParameter(u,deltaR[u],deltaRErr[u]);
    std::cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << std::endl;
    myMin2->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    std::cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << std::endl;
  }

  std::cout << " ################## second step done ##############" << std::endl;

  // TMinuit *myMin3 = new TMinuit(192);
  // myMin3->SetObjectFit(quickOptAllAntStepsVPOL);
  // myMin3->SetFCN(iHateRoot);

  // for(int y = 0; y <MAX_ANTENNAS; y++){    
  //   char name[30];
  //   sprintf(name,"r%d",y);
  //   myMin3->DefineParameter(y, name, deltaR[y], stepSize, -0.3,0.3);
  //   sprintf(name,"z%d",y);
  //   myMin3->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
  //   sprintf(name,"phi%d",y);
  //   myMin3->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
  //   sprintf(name,"cable%d",y);
  //   myMin3->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
  // }



  // for(int y = 0; y <MAX_ANTENNAS; y++){
  //   // myMin3->FixParameter(y); // fixed R 
  //   myMin3->FixParameter(y+MAX_ANTENNAS); // fixed Z
  //   myMin3->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
  //   myMin3->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  // }

  // myMin3->Migrad();


  // for(int u = 0; u <MAX_ANTENNAS; u++){
  //   // myMin3->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
  //   // std::cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << std::endl;
  //   myMin3->GetParameter(u,deltaR[u],deltaRErr[u]);
  //   std::cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << std::endl;

  // }

  // std::cout << " ################## third step done ##############" << std::endl;




  TMinuit *myMin2again = new TMinuit(192);
  myMin2again->SetObjectFit(quickOptAllAntStepsVPOL);
  myMin2again->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin2again->DefineParameter(y, name, deltaR[y], stepSize, -0.3, 0.3);
    sprintf(name,"z%d",y);
    myMin2again->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
    sprintf(name,"phi%d",y);
    myMin2again->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
    sprintf(name,"cable%d",y);
    myMin2again->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
    myMin2again->FixParameter(y); // fixed R 
    myMin2again->FixParameter(y+MAX_ANTENNAS); // fixed Z
    myMin2again->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
    // myMin2again->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  myMin2again->Migrad();


  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin2again->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    std::cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << std::endl;

  }


  std::cout << " ################second half step done #############" << std::endl;



  TMinuit *myMin4 = new TMinuit(192);
  myMin4->SetObjectFit(quickOptAllAntStepsVPOL);
  myMin4->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin4->DefineParameter(y, name, deltaR[y],stepSize, -0.3, 0.3);
    sprintf(name,"z%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
    sprintf(name,"phi%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
    sprintf(name,"cable%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
     myMin4->FixParameter(y); // fixed R 
     //myMin4->FixParameter(y+MAX_ANTENNAS); // fixed Z
     myMin4->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
     myMin4->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin4->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  


  std::time_t now = std::time(NULL);
  std::tm * ptm = std::localtime(&now);
  char buffer[32];
  std::strftime(buffer, 32, "%Y_%m_%d_time_%H_%M_%S", ptm); 


  std::ofstream newfile(Form("RETEST/newLindaNumbers_4steps_VPOL_10kVSeavey_NEW14_onlyMean_%s.txt", buffer));
  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin4->GetParameter(u,deltaR[u],deltaRErr[u]);
    std::cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << std::endl;
	
    
    myMin4->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    std::cout << " deltaZ[" << u << "] = " << deltaZ[u] << " +/- " << deltaZErr[u] << std::endl;
    
    myMin4->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    std::cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << std::endl;

    myMin4->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    std::cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << std::endl;
    
    newfile << u << "  " << deltaR[u]<< "  " << deltaZ[u]<< "  " << deltaPhi[u]<< "  " << deltaCableDelays[u] << std::endl;


  }
  
  
  std::cout << "Easy table" << std::endl;
  for(int u = 0; u <MAX_ANTENNAS; u++)  std::cout << u << " & " << deltaR[u]<< " & " << deltaZ[u]<< " & " << deltaPhi[u]<< " & " << deltaCableDelays[u] << std::endl;

  
}

