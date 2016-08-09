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

void runQuickDeltazVPOL() {
  
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  
  double startVal=0;               
  double stepSize=0.0000005;       
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
  gSystem->CompileMacro("quickDeltazVPOL.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  Double_t relDeltaOut=0;

  const int ndt = 30;

  Double_t arrayDeltaT[ndt];
  Double_t arrayDeltaZ[ndt];

  for (int idelta=0; idelta<ndt;idelta++){
    
    arrayDeltaT[idelta]= -0.2 + (0.6/ndt)*idelta;
  
  TMinuit *myMin = new TMinuit(3);
  myMin->SetObjectFit(quickDeltazVPOL);
  myMin->SetFCN(iHateRoot);
 //   myMin->SetMaxIterations(2);
  //myMin->SetErrorDef(1000);

//   int ierflg;
//   double eps[1] = {2.};
//   myMin->mnexcm("SET EPS", eps, 1, ierflg);
  //setArray();

  double deltaZ=0;
  double deltaZErr;
  double deltaR=0;
  double deltaRErr;
  double deltaPHI=0;
  double deltaPHIErr;
  double deltaT=arrayDeltaT[idelta];
  double deltaTErr;


  char name[30];
  sprintf(name,"deltaZ");
  myMin->DefineParameter(0, name, deltaZ, stepSize, -1, 1);    
  // sprintf(name,"deltaZ");
  // myMin->DefineParameter(1, name, deltaZ, stepSize, minVal, maxVal);    
  // sprintf(name,"deltaPHI");
  // myMin->DefineParameter(2, name, deltaPHI, stepSize, minVal, maxVal);    
  sprintf(name,"deltaT");
   myMin->DefineParameter(1, name, deltaT, stepSize, minVal, maxVal);    

   myMin->FixParameter(1);

  //*********MINUIT METHOD*******************
//   myMin->SetPrintLevel(-1);
  myMin->Migrad();
  //  int error_flag;
  //  myMin->mnexcm("MINOS",0,0,error_flag);  

  myMin->GetParameter(0,deltaZ,deltaZErr);
  cout << " deltaZ[] = " << deltaZ << " +/- " << deltaZErr << endl;
  // myMin->GetParameter(1,deltaZ,deltaZErr);
  // cout << " deltaZ[] = " << deltaZ << " +/- " << deltaZErr << endl;    
  // myMin->GetParameter(2,deltaPHI,deltaPHIErr);
  // cout << " deltaPHI[] = " << deltaPHI << " +/- " << deltaPHIErr << endl;    
  myMin->GetParameter(1,deltaT,deltaTErr);
  cout << " deltaT[] = " << deltaT << " +/- " << deltaTErr << endl; 

  arrayDeltaZ[idelta] = deltaZ;
  }




  for (int i=0;i<ndt;i++){
    cout << arrayDeltaT[i] << " " << arrayDeltaZ[i] << endl;
  }



  TCanvas *c1 = new TCanvas("c1");
  TGraph *g = new TGraph(ndt, arrayDeltaT, arrayDeltaZ);
  g->SetTitle(";Cable delay [ns]; #Delta z [m]");

  g->SetMarkerColor(kBlue);
  g->SetMarkerStyle(4);
  g->Draw("Apl");


  TLine * l1 = new TLine(0, c1->GetUymin(), 0, c1->GetUymax());
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);
  l1->Draw();

  TLine * l2 = new TLine(c1->GetUxmin(), arrayDeltaZ[10], c1->GetUxmax(), arrayDeltaZ[10]);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(2);
  l2->Draw();

  TLine * l3 = new TLine(arrayDeltaT[20], c1->GetUymin(), arrayDeltaT[20], c1->GetUymax());
  l3->SetLineColor(kViolet);
  l3->SetLineStyle(2);
  l3->Draw();

  TLine * l4 = new TLine(c1->GetUxmin(), 0, c1->GetUxmax(), 0);
  l4->SetLineColor(kViolet);
  l4->SetLineStyle(2);
  l4->Draw();

  c1->Print("VPOL_dt_dz_correlation.png");
  c1->Print("VPOL_dt_dz_correlation.pdf");
  

}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){

  double diffErr = quickDeltazVPOL(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


