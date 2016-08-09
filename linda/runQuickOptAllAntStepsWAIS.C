//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48

const int reallyBigNumber = 2300000;
double eventNumberIndex[reallyBigNumber];
double thetaWaveIndex[reallyBigNumber];
double phiWaveIndex[reallyBigNumber];
int antIndex1[reallyBigNumber];
int antIndex2[reallyBigNumber];
double maxCorrTimeIndex[reallyBigNumber];
bool adjacent[reallyBigNumber];

void runQuickOptAllAntStepsWAIS() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.0001;       
  double minVal=-0.5;           
  double maxVal=0.5;            
  Double_t p0 = 0;

  //Load libraries. Need to have ANITA_UTIL_INSTALL_DIR/lib and ROOTSYS/lib in the LD_LIBRARY_PATH               
  // gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");      
  gSystem->Load("libAnitaCorrelator.so");

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->useKurtAnita3Numbers(1);
  gSystem->CompileMacro("quickOptAllAntStepsWAIS.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(192);
  myMin->SetObjectFit(quickOptAllAntStepsWAIS);
  myMin->SetFCN(iHateRoot);
//   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();


  Double_t deltaR[MAX_ANTENNAS],deltaRErr[MAX_ANTENNAS];
  Double_t deltaZ[MAX_ANTENNAS],deltaZErr[MAX_ANTENNAS];
  Double_t deltaPhi[MAX_ANTENNAS],deltaPhiErr[MAX_ANTENNAS];
  Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];


  for(int y = 0; y <MAX_ANTENNAS; y++){

    deltaR[y] = 0.;
    deltaZ[y] = 0.;
    deltaPhi[y] = 0.;
    deltaCableDelays[y] = 0.;
    
    if (y==47) deltaCableDelays[y] = -0.35;

    char name[30];
    sprintf(name,"r%d",y);
    myMin->DefineParameter(y, name, deltaR[y], stepSize, -0.15, 0.15);
    sprintf(name,"z%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.15, 0.15);
    sprintf(name,"phi%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.15, 0.15);
    sprintf(name,"cable%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, deltaCableDelays[y]-0.15, deltaCableDelays[y]+0.15 );
  }
  
  
  
  Double_t deltaR[MAX_ANTENNAS],deltaRErr[MAX_ANTENNAS];
  Double_t deltaZ[MAX_ANTENNAS],deltaZErr[MAX_ANTENNAS];
  Double_t deltaPhi[MAX_ANTENNAS],deltaPhiErr[MAX_ANTENNAS];
  Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];
  

  myMin->FixParameter(0+MAX_ANTENNAS*3); // fixed t0

  for(int y = 0; y <MAX_ANTENNAS; y++){
    myMin->FixParameter(y); // fixed R 
    myMin->FixParameter(y+MAX_ANTENNAS); // fixed Z
    myMin->FixParameter(y+MAX_ANTENNAS*2); // fixed phi

  }

  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  



  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << endl;
  }

  cout << " ################## first step done ##############" << endl;
  
  TMinuit *myMin2 = new TMinuit(192);
  myMin2->SetObjectFit(quickOptAllAntStepsWAIS);
  myMin2->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin2->DefineParameter(y, name, deltaR[y], stepSize, -0.15, 0.15);
    sprintf(name,"z%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.15, 0.15);
    sprintf(name,"phi%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.02, 0.02);
    sprintf(name,"cable%d",y);
    myMin2->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, -0.5, 0.5);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
    myMin2->FixParameter(y); // fixed R 
    myMin2->FixParameter(y+MAX_ANTENNAS); // fixed Z
    // myMin2->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
    myMin2->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  myMin2->Migrad();


  for(int u = 0; u <MAX_ANTENNAS; u++){
    // myMin2->GetParameter(u,deltaR[u],deltaRErr[u]);
    // cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << endl;
    myMin2->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << endl;
  }


  cout << " ################## second step done ##############" << endl;

  TMinuit *myMin3 = new TMinuit(192);
  myMin3->SetObjectFit(quickOptAllAntStepsWAIS);
  myMin3->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin3->DefineParameter(y, name, deltaR[y], stepSize, -0.15,0.15);
    sprintf(name,"z%d",y);
    myMin3->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.15, 0.15);
    sprintf(name,"phi%d",y);
    myMin3->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.15, 0.15);
    sprintf(name,"cable%d",y);
    myMin3->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, -0.5, 0.5);
  }



  for(int y = 0; y <MAX_ANTENNAS; y++){
    //    myMin3->FixParameter(y); // fixed R 
    myMin3->FixParameter(y+MAX_ANTENNAS); // fixed Z
    myMin3->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
    myMin3->FixParameter(y+MAX_ANTENNAS*3); // fixed t

  }

  myMin3->Migrad();


  for(int u = 0; u <MAX_ANTENNAS; u++){
    //   myMin3->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    //   cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << endl;
    myMin3->GetParameter(u,deltaR[u],deltaRErr[u]);
    cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << endl;

  }

  cout << " ################## third step done ##############" << endl;

  TMinuit *myMin4 = new TMinuit(192);
  myMin4->SetObjectFit(quickOptAllAntStepsWAIS);
  myMin4->SetFCN(iHateRoot);

  for(int y = 0; y <MAX_ANTENNAS; y++){    
    char name[30];
    sprintf(name,"r%d",y);
    myMin4->DefineParameter(y, name, deltaR[y],stepSize, -0.15, 0.15);
    sprintf(name,"z%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.15, 0.15);
    sprintf(name,"phi%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.15, 0.15);
    sprintf(name,"cable%d",y);
    myMin4->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, -0.5, 0.5);
  }


  for(int y = 0; y <MAX_ANTENNAS; y++){
    myMin4->FixParameter(y); // fixed R 
    //myMin->FixParameter(y+MAX_ANTENNAS); // fixed Z
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


  ofstream newfile(Form("RETEST/newLindaNumbers_4steps_WAISHPOL_NEW12_cosminV3_nfixedBug_withOffAxisDelay_%s.txt", buffer));
  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin4->GetParameter(u,deltaR[u],deltaRErr[u]);
    cout << "deltaR[" << u << "] = " << deltaR[u] << " +/- " << deltaRErr[u] << endl;
	
    
    myMin4->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    cout << " deltaZ[" << u << "] = " << deltaZ[u] << " +/- " << deltaZErr[u] << endl;
    
    myMin4->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << " +/- " << deltaPhiErr[u] << endl;

    myMin4->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << endl;
    
    newfile << u << "  " << deltaR[u]<< "  " << deltaZ[u]<< "  " << deltaPhi[u]<< "  " << deltaCableDelays[u] << endl;


  }
  
  
  cout << "Easy table" << endl;
  for(int u = 0; u <MAX_ANTENNAS; u++)  cout << u << " & " << deltaR[u]<< " & " << deltaZ[u]<< " & " << deltaPhi[u]<< " & " << deltaCableDelays[u] << endl;

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = quickOptAllAntStepsWAIS(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


