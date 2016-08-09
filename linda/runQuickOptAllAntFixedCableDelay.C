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

void runQuickOptAllAntFixedCableDelay() {
  
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
  gSystem->CompileMacro("quickOptAllAntFixedCableDelay.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(192);
  myMin->SetObjectFit(quickOptAllAntFixedCableDelay);
  myMin->SetFCN(iHateRoot);
//   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

  double startValR[MAX_ANTENNAS] ={0};
  double startValZ[MAX_ANTENNAS] ={0};
  double startValPhi[MAX_ANTENNAS] ={0};
  double startValCableDelays[MAX_ANTENNAS] ={0};

  int tempint;
  ifstream input("newLindaNumbersCableDelaysONLY.txt");
  for(int y = 0; y <MAX_ANTENNAS; y++){
    input >> tempint >> startValCableDelays[y];
    //   cout << tempint << " " << startValCableDelays[y] << endl;
  }

  input.close();
  
  for(int y = 0; y <MAX_ANTENNAS; y++){
    
    char name[30];
    sprintf(name,"r%d",y);
    myMin->DefineParameter(y, name, startValR[y], stepSize, -0.15, 0.001);
    sprintf(name,"z%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS, name, startValZ[y], stepSize, -0.1, 0.1);
    sprintf(name,"phi%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*2, name, startValPhi[y], stepSize, -0.1, 0.1);
    sprintf(name,"cable%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*3, name, startValCableDelays[y], stepSize, -0.2, 0.2);
    myMin->FixParameter(y+MAX_ANTENNAS*3);
  }
  
  
  
  Double_t deltaR[MAX_ANTENNAS],deltaRErr[MAX_ANTENNAS];
  Double_t deltaZ[MAX_ANTENNAS],deltaZErr[MAX_ANTENNAS];
  Double_t deltaPhi[MAX_ANTENNAS],deltaPhiErr[MAX_ANTENNAS];
  Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];
  
  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  

  ofstream newfile("newLindaNumbersCableDelaysFixed_ALLANTENNAS.txt");
  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin->GetParameter(u,deltaR[u],deltaRErr[u]);
    cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << endl;
	
    
    myMin->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    cout << " deltaZ[" << u << "] = " << deltaZ[u] << " +/- " << deltaZErr[u] << endl;
    
    myMin->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << endl;

    myMin->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << endl;
    
    newfile << u << "  " << deltaR[u]<< "  " << deltaZ[u]<< "  " << deltaPhi[u]<< "  " << deltaCableDelays[u] << endl;


  }
  
  
  cout << "Easy table" << endl;
  for(int u = 0; u <MAX_ANTENNAS; u++)  cout << u << " & " << deltaR[u]<< " & " << deltaZ[u]<< " & " << deltaPhi[u]<< " & " << deltaCableDelays[u] << endl;

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = quickOptAllAntFixedCableDelay(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


