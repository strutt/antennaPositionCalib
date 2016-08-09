//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48

void runTiltAllAnt() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.1;       
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
  gSystem->CompileMacro("tiltAllAnt.C","k");
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(2);
  myMin->SetObjectFit(tiltAllAnt);
  myMin->SetFCN(iHateRoot);
  //myMin->SetErrorDef(1000);

  double tilt[2];
  double tiltErr[2];

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

//   double startValDeltaT[MAX_ANTENNAS] ={0};
//   double startValR[MAX_ANTENNAS] ={0};
//   double startValPhi[MAX_ANTENNAS] ={0};
//   double startValCableDelays[MAX_ANTENNAS] ={0};
  
  
//   for(int y = 0; y <MAX_ANTENNAS; y++){
    
//     char name[30];
//     sprintf(name,"r%d",y);
    myMin->DefineParameter(0, "tilt", 0, 0.01, -0.1, 0.1);
    myMin->DefineParameter(0, "tiltPos", 0, 0.1, -TMath::Pi(), TMath::Pi());
//     sprintf(name,"z%d",y);
//     myMin->DefineParameter(y+MAX_ANTENNAS, name, startValDeltaT[y], stepSize, minVal, maxVal);
//     sprintf(name,"phi%d",y);
//     myMin->DefineParameter(y+MAX_ANTENNAS*2, name, startValPhi[y], stepSize, minVal, maxVal);
//     sprintf(name,"cable%d",y);
//     myMin->DefineParameter(y+MAX_ANTENNAS*3, name, startValCableDelays[y], stepSize, minVal, maxVal);
//   }
  
  
  
//   Double_t deltaR[MAX_ANTENNAS],deltaRErr[MAX_ANTENNAS];
//   Double_t deltaZ[MAX_ANTENNAS],deltaZErr[MAX_ANTENNAS];
//   Double_t deltaPhi[MAX_ANTENNAS],deltaPhiErr[MAX_ANTENNAS];
//   Double_t deltaCableDelays[MAX_ANTENNAS],deltaCableDelaysErr[MAX_ANTENNAS];
  
  //*********MINUIT METHOD*******************
  myMin->SetPrintLevel(-1);
  myMin->Migrad();
  

  myMin->GetParameter(0,tilt[0],tiltErr[0]);
  myMin->GetParameter(1,tilt[1],tiltErr[1]);


  cout << "TILT : " << tilt[0] << " " << tiltErr[0] << endl;
  cout << "TILT POS : " << tilt[1] << " " << tiltErr[1] << endl;


  
  
  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = tiltAllAnt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


