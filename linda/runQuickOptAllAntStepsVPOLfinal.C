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

void runQuickOptAllAntStepsVPOLfinal() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.0001;       
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

  gSystem->CompileMacro("quickOptAllAntStepsVPOL.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

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


  // ifstream myInput("RETEST/newLindaNumbers_4steps_VPOL_10kVSeavey_NEW12_bugFixed_2016_02_12_time_10_47_14.txt");
  ifstream myInput("RETEST/newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt");

  int tempCount = 0;
  int random;
   if (myInput.is_open()){
     while(tempCount <48){
       myInput >> random >> deltaR[tempCount] >> deltaZ[tempCount] >> deltaPhi[tempCount] >> deltaCableDelays[tempCount];

       cout <<  random << " " << deltaR[tempCount] << " " <<  deltaZ[tempCount] << " " <<  deltaPhi[tempCount] << " "<< deltaCableDelays[tempCount] << endl;
       // deltaZ[tempCount] =0;
       tempCount++;
     }   
   }
  
  for(int y = 0; y <MAX_ANTENNAS; y++){
 
    minValCableDelays[y] = -0.5;
    maxValCableDelays[y] = +0.5;
    
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
    cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u]  << " +/- " << deltaCableDelaysErr[u] << endl;

    //    myMin->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    //    cout << " deltaZ[" << u << "] = " << deltaZ[u] << " +/- " << deltaZErr[u] << endl;

  }

  cout << " ################## third step done ##############" << endl;

//   TMinuit *myMin4 = new TMinuit(192);
//   myMin4->SetObjectFit(quickOptAllAntStepsVPOL);
//   myMin4->SetFCN(iHateRoot);

//   for(int y = 0; y <MAX_ANTENNAS; y++){    

//     sprintf(name,"r%d",y);
//     myMin4->DefineParameter(y, name, deltaR[y],stepSize, -0.3, 0.3);
//     sprintf(name,"z%d",y);
//     myMin4->DefineParameter(y+MAX_ANTENNAS, name, deltaZ[y], stepSize, -0.3, 0.3);
//     sprintf(name,"phi%d",y);
//     myMin4->DefineParameter(y+MAX_ANTENNAS*2, name, deltaPhi[y], stepSize, -0.3, 0.3);
//     sprintf(name,"cable%d",y);
//     myMin4->DefineParameter(y+MAX_ANTENNAS*3, name, deltaCableDelays[y], stepSize, minValCableDelays[y], maxValCableDelays[y]);
//   }


//   for(int y = 0; y <MAX_ANTENNAS; y++){
//      myMin4->FixParameter(y); // fixed R 
//      //myMin4->FixParameter(y+MAX_ANTENNAS); // fixed Z
//      myMin4->FixParameter(y+MAX_ANTENNAS*2); // fixed phi
//      myMin4->FixParameter(y+MAX_ANTENNAS*3); // fixed t

//   }

//   //*********MINUIT METHOD*******************
// //   myMin->SetPrintLevel(-1);
//   myMin4->Migrad();
//   //  int error_flag;
//   //  myMin->mnexcm("MINOS",0,0,error_flag);  


  std::time_t now = std::time(NULL);
  std::tm * ptm = std::localtime(&now);
  char buffer[32];
  std::strftime(buffer, 32, "%Y_%m_%d_time_%H_%M_%S", ptm); 


  ofstream newfile(Form("RETEST/newLindaNumbers_4steps_VPOL_10kVSeavey_NEW12_bugFixed_reFitFromWAISGEOMETRY_%s.txt", buffer));
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

  double diffErr = quickOptAllAntStepsVPOL(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


