//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48

void runOptAllAnt() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.05;       
  double minVal=-0.1;           
  double maxVal=0.1;            
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
  gSystem->CompileMacro("optAllAnt.C","k");
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(192);
  myMin->SetObjectFit(optAllAnt);
  myMin->SetFCN(iHateRoot);
//   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

  double startValDeltaT[MAX_ANTENNAS] ={0};
  double startValR[MAX_ANTENNAS] ={0};
  double startValPhi[MAX_ANTENNAS] ={0};
  double startValCableDelays[MAX_ANTENNAS] ={0};
  
  
  for(int y = 0; y <MAX_ANTENNAS; y++){
    
    char name[30];
    sprintf(name,"r%d",y);
    myMin->DefineParameter(y, name, startValR[y], stepSize, minVal, maxVal);
    sprintf(name,"z%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS, name, startValDeltaT[y], stepSize, minVal, maxVal);
    sprintf(name,"phi%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*2, name, startValPhi[y], stepSize, minVal, maxVal);
    sprintf(name,"cable%d",y);
    myMin->DefineParameter(y+MAX_ANTENNAS*3, name, startValCableDelays[y], stepSize, minVal, maxVal);
  }
  
  
  
  Double_t deltaR[MAX_ANTENNAS],deltaRErr[MAX_ANTENNAS];
  Double_t deltaZ[MAX_ANTENNAS],deltaZErr[MAX_ANTENNAS];
  Double_t deltaPhi[MAX_ANTENNAS],deltaPhiErr[MAX_ANTENNAS];
  Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];
  
  //*********MINUIT METHOD*******************
  myMin->SetPrintLevel(-1);
  myMin->Migrad();
//   int error_flag;
//   myMin->mnexcm("SIMPLEX",0,0,error_flag);  

  ofstream newfile("newLindaNumbersCableDelays_ALLANTENNAS_GRAD.txt");
  for(int u = 0; u <MAX_ANTENNAS; u++){
    myMin->GetParameter(u,deltaR[u],deltaRErr[u]);
    //cout << "deltaR[" << u << "] = " << deltaR[u] ;
	
    
    myMin->GetParameter(u+MAX_ANTENNAS,deltaZ[u],deltaZErr[u]);
    //cout << " deltaZ[" << u << "] = " << deltaZ[u] ;
    
    myMin->GetParameter(u+MAX_ANTENNAS*2,deltaPhi[u],deltaPhiErr[u]);
    //cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << ";" << endl;

    myMin->GetParameter(u+MAX_ANTENNAS*3,deltaCableDelays[u],deltaCableDelaysErr[u]);
    //cout << " deltaCableDelays[" << u << "] = " << deltaCableDelays[u] << ";" << endl;
    
    newfile << u << "  " << deltaZ[u]<< "  " << deltaR[u]<< "  " << deltaPhi[u]<< "  " << deltaCableDelays[u] << endl;


  }
  
  
  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = optAllAnt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


