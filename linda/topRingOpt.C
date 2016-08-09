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


double topRingOpt(double *par);
vector<double> leastSquares(double *xIn, double *yIn, int count);
double getMean(double *x, int n);
double getRMS(double *x, int n);

double topRingOpt(double *par){


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

  Double_t deltaR[48]={0};  
  Double_t deltaZ[48]={0};  
  Double_t deltaPhi[48]={0};  
  Double_t deltaHeading[1]={0};
  double  deltaTArrayMod[48]={0};
 

  for(int i = 0; i<16;i++){
    deltaR[i]=par[i];
    deltaZ[i]=par[i+16];
    deltaPhi[i]=par[i+32];
    //deltaTArrayMod[i] = par[i+16];
  }

 
  double theReturn = 0;
  double sumMean = 0;
  double sumRMS = 0;
  double sumMean2 = 0;
  int count8 = 0;
  double sumGrads = 0;


  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  char *corTreeDir = "../corTrees/";
  double dummyArray[MAX_ANTENNAS][1] ={{0}}; 
  TGraph *tempAntGraph;

  vector<vector<double> > phiAngle;
  vector<vector<double> > deltaTVec;
  vector<vector<int> > firstAntVec;
  vector<vector<int> > secondAntVec;

  vector<vector<double> > phiAngleArray2;
  vector<vector<double> > deltaTArray2;

  vector<double> temp;
  vector<int> temp2;
  temp.push_back(0);
  temp2.push_back(0);

  double deltaTArrayLoop[6000] ={0};
  double phiAngleArrayLoop[6000] = {0};


  int leftOpt, rightOpt;

  double meanPhi[48] = {0}; 
  for (int ant=0;ant<48;ant++){
    meanPhi[ant] = fGeomTool->getAntPhiPositionRelToAftFore(ant, pol);
    //    meanPhi[ant] = fGeomTool->getAntPhiPosition(ant, pol);
  }

  Double_t additionalPhi = 22.5*TMath::DegToRad();

  for(int i =0; i < MAX_ANTENNAS; i++){
    phiAngleArray2.push_back(temp);
    deltaTArray2.push_back(temp);

  }

  for(int run = 332; run <356; run++){ // only 1 run for now
    

    phiAngle.clear();
    deltaTVec.clear();
    firstAntVec.clear();
    secondAntVec.clear();

    for(int i = 0; i < MAX_ANTENNAS; i++){
      phiAngle.push_back(temp);
      deltaTVec.push_back(temp);
      firstAntVec.push_back(temp2);
      secondAntVec.push_back(temp2);
    }

    sprintf(baseDir,"/unix/anita3/flight1415/root");
    sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
    sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
    sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
    sprintf(corrName,"%s/corRun_NEW2_HPOL_%d.root",corTreeDir,run);

   
    RawAnitaEvent *event = 0;
    PrettyAnitaHk *hk = 0;
   
    RawAnitaHeader *header =0;
    Adu5Pat *pat =0;
    CorrelationSummaryAnita3 *corSum =0;
   
  
    TFile *fpHead = TFile::Open(headerName);
    TTree *headTree = (TTree*) fpHead->Get("headTree");
    headTree->SetBranchAddress("header",&header);
    headTree->BuildIndex("eventNumber");
   
    TFile *fpGps = TFile::Open(gpsName);
    TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
    adu5PatTree->BuildIndex("realTime");
    adu5PatTree->SetBranchAddress("pat",&pat);
   
//     Int_t labChip;
    TFile *fpCor = new TFile(corrName);
    TTree *corTree = (TTree*) fpCor->Get("corTree");
    corTree->SetBranchAddress("cor",&corSum);
//     corTree->SetBranchAddress("labChip",&labChip);

    Long64_t numEntries=corTree->GetEntries();
    int counter=0;

    Long64_t entry=0;
    UInt_t eventNumber, triggerTime, triggerTimeNs;
    Int_t firstAnt,secondAnt,maxAnt,corInd;
    Double_t deltaT,deltaTExpected;
    Double_t phiWave, phiMaxAnt;
    Double_t corPeak, corRMS;
    Double_t balloonLat, balloonLon, balloonAlt;
    Double_t heading,pitch,roll;
  

    Double_t thetaWave;

    for(entry=0;entry<numEntries;entry++) {

      corTree->GetEntry(entry);
      Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
      if(headEntry<0) continue;
      headTree->GetEntry(headEntry);
     
      triggerTimeNs=header->triggerTimeNs;
      triggerTime=header->triggerTime;
      //      eventNumber=header->eventNumber;
  

      Long64_t bestEntry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
      if(bestEntry>-1) adu5PatTree->GetEntry(bestEntry);
      else continue;
      
      balloonLat=pat->latitude;
      balloonLon=pat->longitude;
      balloonAlt=pat->altitude;
      heading=pat->heading;
//       pat->pitch=0.64;
//       pat->roll=0.14;      
//       pitch=pat->pitch;
//       roll=pat->roll;

      UsefulAdu5Pat usefulPat(pat);
      
     for(int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;corInd++) {

       if (corSum->maxCorVals[corInd]/corSum->rmsCorVals[corInd]<6) continue;

       if (corInd<5 || corInd>12) continue; // only neighbouring antennas
       
       deltaTExpected=usefulPat.getDeltaTExpectedOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd], sourceLon, sourceLat, sourceAlt, deltaR, deltaZ, deltaPhi);

	
	firstAnt=corSum->firstAnt[corInd];
	secondAnt=corSum->secondAnt[corInd];

	deltaT=corSum->maxCorTimes[corInd];
	maxAnt=corSum->centreAntenna;
	phiWave=usefulPat.getPhiWave();
	phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	corPeak=corSum->maxCorVals[corInd];
	corRMS=corSum->rmsCorVals[corInd];


	if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 1 ){
	  phiAngle[0].push_back(phiWave);
	  deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
	  firstAntVec[0].push_back(firstAnt);
	  secondAntVec[0].push_back(secondAnt);  
	}
      }

      counter++; 
    
    }

 

    double deltaTArray[48][6000] = {{0}};
    double phiAngleArray[48][6000]= {{0}};

    double deltaTArrayCut[48][6000]= {{0}};
    double phiAngleArrayCut[48][6000]= {{0}};

    int middleAnt; 
    int leftAnt,rightAnt;

    int countArray[48] = {0};


    //fill arrays
    //for(int ants = par[0]; ants < par[0]+1; ants++){
    for(int ants = 0; ants < 16; ants++){

      int count = 0;
      int count2 = 0;
      double sumPhi = 0;
      bool true1 = false;
      bool true2 = false;

      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      
      for(int events = 1; events < phiAngle[0].size(); events++){
	int firstAntTemp = (int)firstAntVec[0][events];
	int secondAntTemp = (int)secondAntVec[0][events];
	int rightTemp = int(rightAnt);

	
	if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	  deltaTArray[ants][count] = deltaTVec[0][events];
	  phiAngleArray[ants][count] = phiAngle[0][events]; 
	  count++;
	}
	
    
      }
      countArray[ants] = count;
      
    }
    
    
    //make cuts
    for(int ants = 0; ants < 16; ants++){
      int count = 0;
      
      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      

      double sumPhi = 0;

      double lower  = meanPhi[ants] - additionalPhi ;
      double upper = meanPhi[rightAnt] + additionalPhi;
      if (lower<0) lower+=2*TMath::Pi();
      if (upper>2*TMath::Pi()) upper-=2*TMath::Pi();
      
      //      cout << lower << " " << upper << endl;


	for(int events = 0; events < countArray[ants]; events++){

	  if (lower>upper){
	    if (phiAngleArray[ants][events]<TMath::Pi()) lower-=2*TMath::Pi();
	    else upper+=2*TMath::Pi();
	  }
	  
	  if((phiAngleArray[ants][events] > lower ) && (phiAngleArray[ants][events]< upper)){
	    phiAngleArrayCut[ants][count] = phiAngleArray[ants][events];
	    deltaTArrayCut[ants][count] = deltaTArray[ants][events];
	    
	    count++;
	      
	  }
	}
	
	for(int events = 0; events < count-1; events++){
	phiAngleArray2[ants].push_back(phiAngleArrayCut[ants][events]);
	deltaTArray2[ants].push_back(deltaTArrayCut[ants][events]);
	}
	
    }
    
    delete event;
    delete hk; 
    delete header;
    delete pat;
    delete corSum;

    delete fpHead;
    delete fpGps ;
    delete fpCor;

  }

  sumMean = 0;
  sumRMS = 0;
  sumMean2 = 0;
  sumGrads = 0;
  for(int ants = 0; ants < 16; ants++){
    count8 = 0;
    

    for(int events = 1; events < phiAngleArray2[ants].size(); events++){

      if( deltaTArrayLoop[count8]<1){      
	deltaTArrayLoop[count8] = deltaTArray2[ants][events];
	phiAngleArrayLoop[count8] = phiAngleArray2[ants][events];
	count8++;
      }
      
    }  
    
      if(count8<3) continue;
      
      if(ants <16){  
	
	double mean = getMean(deltaTArrayLoop, count8);
	sumMean += mean*mean;

	double rms = getRMS(deltaTArrayLoop, count8);
	sumRMS += rms;

//   	vector<double> myFit = leastSquares(phiAngleArrayLoop, deltaTArrayLoop, count8-1);
//  	sumGrads += myFit[0]*myFit[0];
      
      }
    
    
  }
    


  
  cout << sumMean << "  " << sumRMS <<endl;
  
  theReturn = sumMean+sumRMS;

  return theReturn;
   
}

vector<double> leastSquares(double *xIn, double *yIn, int count) {
  
  double SUMx, SUMy, SUMxy, SUMxx, SUMres, res, slope,
    y_intercept, y_estimate ;
  int i,n;
  
  SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
  for (i=0; i<count; i++) {
    
    SUMx = SUMx + xIn[i];
    SUMy = SUMy + yIn[i];
    SUMxy = SUMxy + xIn[i]*yIn[i];
    SUMxx = SUMxx + xIn[i]*xIn[i];
  }
  
  slope = ( SUMx*SUMy - count*SUMxy ) / ( SUMx*SUMx - count*SUMxx );
  y_intercept = ( SUMy - slope*SUMx ) / count;


  vector<double> theReturn2;
  theReturn2.push_back(slope);
  theReturn2.push_back(y_intercept);
  return theReturn2;
  
}




double getMean(double *x, int n){

  double mean = 0;
  for (int i=0;i<n;i++) mean+=x[i];

  return mean/n;

}

double getRMS(double *x, int n){ 

  double rms = 0;
  for (int i=0;i<n;i++) rms+=x[i]*x[i];

  return rms/n;

}
