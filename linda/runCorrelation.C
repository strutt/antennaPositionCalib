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

void runCorrelation() {
  
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
  gSystem->CompileMacro("findCorrelation.C","k");

  fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  
  double startValZ[MAX_ANTENNAS] ={0};
  double startValCableDelays[MAX_ANTENNAS] ={0};
  
  int npar = 20;

  double limitz = 0.1;
  double limitt = 0.1;



  double par1[MAX_ANTENNAS*2] = {0};
  double baseline = findCorrelation(par1, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);
  double minFun = 0;


  for (int antenna =32; antenna<MAX_ANTENNAS; antenna++){
    //  int antenna = 2;
  TH2D *hCorr = new TH2D ("hCorr", "", npar, -limitz, limitz, npar, -limitt, limitt);
  cout << "ant : " << antenna << endl;
  double par[MAX_ANTENNAS*2] = {0};

  for (int i=0; i<npar; ++i){ // index on z
    for (int j=0; j<npar; ++j){ // index on cable delays
    
      par[antenna] = -limitz + (2*limitz)*(i+0.5)*1./npar;
      par[antenna+MAX_ANTENNAS]= -limitt + (2*limitt)*(j+0.5)*1./npar;

      minFun = findCorrelation(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent);

      cout << par[antenna] << " " << par[antenna+MAX_ANTENNAS] << " " << baseline-minFun << endl;

      hCorr->Fill(par[antenna], par[antenna+MAX_ANTENNAS], minFun-baseline);

    }
  }



  const int NRGBs = 3, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
  Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
  Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
  Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1");

  hCorr->SetTitle(Form("Antenna %i;#Delta z [m];Time offset [ns]", antenna));

  double max = hCorr->GetMaximum();
  if (max<-hCorr->GetMinimum()) max = -hCorr->GetMinimum();

  max = 0.1;

  hCorr->SetMinimum(-max);
  hCorr->SetMaximum(max);


  hCorr->Draw("colz");

  c1->Print(Form("Correlation_dzdt_Ant%i.png", antenna));
  c1->Print(Form("Correlation_dzdt_Ant%i.pdf", antenna));

  }

}

