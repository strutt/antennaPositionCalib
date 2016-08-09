//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

#define MAX_ANTENNAS 48

void runTopRingOpt2() {
  
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
  gSystem->CompileMacro("topRingOpt2.C","k");
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(200);
  myMin->SetObjectFit(topRingOpt2);
  myMin->SetFCN(iHateRoot);
  //setArray();

  double startValDeltaT[MAX_ANTENNAS] ={0};
  double startValR[MAX_ANTENNAS] ={0};
  double startValPhi[MAX_ANTENNAS] ={0};
  double startValCableDelays[MAX_ANTENNAS] ={0};
  
  //    startValDeltaT[0] =  -0.0519515;     startValR[0] =   -0.0101463;         startValPhi[0] =      -0.00473836 ;   
  //     startValDeltaT[1] =     -0.0597062;   startValR[1] =      -0.02577;       startValPhi[1] =   0.00864501 ; 
  //     startValDeltaT[2] =     -0.081435;     startValR[2] =     -0.000224044;    startValPhi[2] =   -0.000630649;
  //     startValDeltaT[3] =     0.0118873;     startValR[3] =      0.019945;       startValPhi[3] =    0.014016;   
  //     startValDeltaT[4] =      0.017917;      startValR[4] =     -0.00297559;    startValPhi[4] =    0.0224936 ;
  //     startValDeltaT[5] =      0.0377119;     startValR[5] =     -0.014872;       startValPhi[5] =   0.0163349; 
  //     startValDeltaT[6] =     -0.0426158;     startValR[6] =    -0.0562555;      startValPhi[6] =    0.0220065 ; 
  //     startValDeltaT[7] =     -0.0221673;     startValR[7] =    -0.034104 ;      startValPhi[7] =    0.0158545 ; 
  //     startValDeltaT[8] =      0.0263739;     startValR[8] =     0.00248804;      startValPhi[8] =    0.013246 ; 
  //     startValDeltaT[9] =      -0.0938419;    startValR[9] =     -0.00344703;     startValPhi[9] =    -0.00718616; 
  //     startValDeltaT[10] =      0.145264;      startValR[10] =     -0.0121874 ;     startValPhi[10] =     0.0156988 ;  
  //     startValDeltaT[11] =      0.118105;      startValR[11] =     -0.0337033 ;      startValPhi[11] =   -0.00324182 ;
  //     startValDeltaT[12] =      0.321805;      startValR[12] =     0.0134362  ;      startValPhi[12] =   -0.00190277 ;
  //     startValDeltaT[13] =      0.0197693;     startValR[13] =     -0.000656063;      startValPhi[13] =   -0.0162318  ;
  //     startValDeltaT[14] =      -0.115263;     startValR[14] =     0.0495637   ;      startValPhi[14] =   -0.0198119 ; 
  //     startValDeltaT[15] =      -0.255707;   startValR[15] =     0.00189892  ;     startValPhi[15] =     0.0383932  ;
  
  
  
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
  //	myMin->SetPrintLevel(-1);
  myMin->Migrad();   
  

  ofstream newfile("newLindaNumbersCableDelays_ALLANTENNAS.txt");
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


  double diffErr = topRingOpt2(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


