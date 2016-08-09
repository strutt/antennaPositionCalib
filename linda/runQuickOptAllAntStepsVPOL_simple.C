//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48


double eventNumberIndex[1500000];
double thetaWaveIndex[1500000];
double phiWaveIndex[1500000];
int antIndex1[1500000];
int antIndex2[1500000];
double maxCorrTimeIndex[1500000];
bool adjacent[1500000];

void runQuickOptAllAntStepsVPOL_simple() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.05;       
  double minVal=-0.5;           
  double maxVal=0.5;            
  Double_t p0 = 0;

  //Load libraries. Need to have ANITA_UTIL_INSTALL_DIR/lib and ROOTSYS/lib in the LD_LIBRARY_PATH               
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");      
  gSystem->Load("libAnitaCorrelator.so");

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->useKurtAnita3Numbers(1);

  gSystem->CompileMacro("quickOptAllAntStepsVPOL_simple.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(192);
  myMin->SetObjectFit(quickOptAllAntStepsVPOL_simple);
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
  

    char name[30];
    // sprintf(name,"r");
    // myMin->DefineParameter(y, name, deltaR[y], stepSize, -0.3, 0.3);
    // sprintf(name,"z");
    // myMin->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
    // sprintf(name,"phi");
    // myMin->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
    sprintf(name,"cable1");
    myMin->DefineParameter(0, name, 0, 0.00001, -0.5, 0.5);
    sprintf(name,"z");
    myMin->DefineParameter(1, name, 0, 0.00001, -0.01, 0.01);
    sprintf(name,"r");
    myMin->DefineParameter(2, name, 0, 0.00001, -0.2, 0.2);



  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  


  
  myMin->GetParameter(0,deltaCableDelays[0],deltaCableDelaysErr[0]);
  
  cout << " deltaCableDelays = " << deltaCableDelays[0]  << " +/- " << deltaCableDelaysErr[0] << endl;

  myMin->GetParameter(1,deltaCableDelays[1],deltaCableDelaysErr[1]);
  
  cout << " deltaZ = " << deltaCableDelays[1]  << " +/- " << deltaCableDelaysErr[1] << endl;

  myMin->GetParameter(2,deltaCableDelays[2],deltaCableDelaysErr[2]);
  
  cout << " deltaCableR = " << deltaCableDelays[2]  << " +/- " << deltaCableDelaysErr[2] << endl;


}


void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = quickOptAllAntStepsVPOL_simple(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


