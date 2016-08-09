//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48


double eventNumberIndex[1300000];
double thetaWaveIndex[1300000];
double phiWaveIndex[1300000];
int antIndex1[1300000];
int antIndex2[1300000];
double maxCorrTimeIndex[1300000];
bool adjacent[1300000];

void runCableDelayOptAllAnt() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.1;       
  double minVal=-0.2;           
  double maxVal=0.2;            
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
  gSystem->CompileMacro("cableOptAllAnt.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(50);
  myMin->SetObjectFit(cableOptAllAnt);
  myMin->SetFCN(iHateRoot);
//   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

  double startValCableDelays[MAX_ANTENNAS] ={0};
  
  
  for(int y = 0; y <MAX_ANTENNAS-1; y++){
    
    char name[30];
    sprintf(name,"cable%d",y);
    myMin->DefineParameter(y, name, startValCableDelays[y], stepSize, minVal, maxVal);
  }

    char name[30];
    sprintf(name,"cable%d",47);
    myMin->DefineParameter(47, name, startValCableDelays[47], stepSize, -0.5, 0.5);  
 
  Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];
  
  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  

  ofstream newfile("newLindaNumbersCableDelaysONLY_2015Oct12.txt");
  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin->GetParameter(u,deltaCableDelays[u],deltaCableDelaysErr[u]);
    cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << endl;
    
    newfile << u << "  " << deltaCableDelays[u] << endl;


  }
  
  
  cout << "Easy table" << endl;
  for(int u = 0; u <MAX_ANTENNAS; u++)  cout << u << " " << deltaCableDelays[u] << endl;

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = cableOptAllAnt(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


