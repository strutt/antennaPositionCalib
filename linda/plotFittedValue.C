#include "AnitaConventions.h"
#include "AnitaGeomTool.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "SurfHk.h"
#include "TurfRate.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>

#define MAX_ANTENNAS 48

void plotFittedValue(){


  AnitaGeomTool *myGeomTool = AnitaGeomTool::Instance();

  double x[MAX_ANTENNAS] = {0.};
  double y[MAX_ANTENNAS] = {0.};
  double z[MAX_ANTENNAS] = {0.};
  double r[MAX_ANTENNAS] = {0.};
  double phi[MAX_ANTENNAS] = {0.};

  double zFit[MAX_ANTENNAS] = {0.};
  double rFit[MAX_ANTENNAS] = {0.};
  double phiFit[MAX_ANTENNAS] = {0.};

  double zFit2[MAX_ANTENNAS] = {0.};
  double rFit2[MAX_ANTENNAS] = {0.};
  double phiFit2[MAX_ANTENNAS] = {0.};


  AnitaPol::AnitaPol_t pol=AnitaPol::kHorizontal;

  //  ifstream myInput("/home/lindac/ANITA/Software/EventCorrelator/macros/lindaMacros/OffsetsErrors_RMSGRAD2_NEWLIM2.txt");
  ifstream myInput("/home/lindac/ANITA/Software/EventCorrelator/macros/lindaMacros/OffsetsErrors_4steps.txt");

   int random;
   int tempCount = 0;
   double temp;

  Double_t deltaR[MAX_ANTENNAS]={0};  
  Double_t deltaZ[MAX_ANTENNAS]={0};  
  Double_t deltaPhi[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays[MAX_ANTENNAS]={0};  

  Double_t deltaRErr[MAX_ANTENNAS]={0};  
  Double_t deltaZErr[MAX_ANTENNAS]={0};  
  Double_t deltaPhiErr[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelaysErr[MAX_ANTENNAS]={0};  

  Double_t deltaR2[MAX_ANTENNAS]={0};  
  Double_t deltaZ2[MAX_ANTENNAS]={0};  
  Double_t deltaPhi2[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays2[MAX_ANTENNAS]={0};  

  Double_t deltaRErr2[MAX_ANTENNAS]={0};  
  Double_t deltaZErr2[MAX_ANTENNAS]={0};  
  Double_t deltaPhiErr2[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelaysErr2[MAX_ANTENNAS]={0};  

  Double_t antArray[MAX_ANTENNAS] = {0.};
  Double_t antErr[MAX_ANTENNAS] = {0.};
  // char *tempChar;

  double tempDouble;


  if (myInput.is_open()){

  for (int i=0;i<48;i++){
    char *tempChar2;
    myInput >> random;
    myInput >> deltaR[i];
    myInput >> deltaRErr[i];
  }
  for (int i=0;i<48;i++) myInput >> random >>  deltaZ[i] >> deltaZErr[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) myInput >> random >>  deltaPhi[i] >> deltaPhiErr[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) myInput >> random >>  deltaCableDelays[i] >> deltaCableDelaysErr[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) cout << i << " " << deltaZ[i] << " " << deltaZErr[i] << endl;
  
  }

  ifstream myInput2("/home/lindac/ANITA/Software/EventCorrelator/macros/lindaMacros/OffsetsErrors_4steps_VPOL.txt");


  if (myInput2.is_open()){

  for (int i=0;i<48;i++){
    char *tempChar2;
    myInput2 >> random;
    myInput2 >> deltaR2[i];
    myInput2 >> deltaRErr2[i];
  }
  for (int i=0;i<48;i++) myInput2 >> random >>  deltaZ2[i] >> deltaZErr2[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) myInput2 >> random >>  deltaPhi2[i] >> deltaPhiErr2[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) myInput2 >> random >>  deltaCableDelays2[i] >> deltaCableDelaysErr2[i];// >> tempDouble >> tempDouble;
  for (int i=0;i<48;i++) cout << i << " " << deltaZ2[i] << " " << deltaZErr2[i] << endl;
  
  }


   for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
     
    myGeomTool->getAntXYZ(ant, x[ant], y[ant], z[ant], pol);
    
    r[ant] = myGeomTool->getAntR(ant, pol);
    phi[ant] = myGeomTool->getAntPhiPositionRelToAftFore(ant, pol);

    cout << ant << "\t" << r[ant] << "\t" << z[ant] << "\t" << phi[ant]*TMath::RadToDeg() << endl;

    rFit[ant] = r[ant] + deltaR[ant];
    zFit[ant] = z[ant] + deltaZ[ant];
    phiFit[ant] = phi[ant] + deltaPhi[ant];
    antArray[ant] = (ant+1)*1.;

    rFit2[ant] = r[ant] + deltaR2[ant];
    zFit2[ant] = z[ant] + deltaZ2[ant];
    phiFit2[ant] = phi[ant] + deltaPhi2[ant];

   }


  string sring[3] = {"Top", "Middle", "Bottom"};

  TCanvas *c1 = new TCanvas("c1", "", 1200, 750);

  TLegend *leg = new TLegend(0.6, 0.11, 0.89, 0.25);
  leg->SetFillColor(kWhite);

 
  TGraph *gphotoR = new TGraph(48, antArray, r);
  gphotoR->SetMarkerStyle(22);
  
  TGraphErrors *gfitR = new TGraphErrors(48, antArray, rFit, antErr, deltaRErr);
  gfitR->SetMarkerStyle(22);
  gfitR->SetMarkerColor(kRed);

  TGraphErrors *gfitR2 = new TGraphErrors(48, antArray, rFit2, antErr, deltaRErr2);
  gfitR2->SetMarkerStyle(22);
  gfitR2->SetMarkerColor(kBlue);

  gphotoR->SetTitle(";Antenna;r [m]");
  
  gphotoR->Draw("Ap");
  gphotoR->SetMinimum(0.);
  gphotoR->SetMaximum(2.6);
  gphotoR->Draw("Ap");

  gfitR->SetLineColor(2);
  gfitR->SetLineWidth(2);
  gfitR->SetFillColor(2);
  gfitR->SetFillStyle(3001);
  gfitR->Draw("p same");

  gfitR2->SetLineColor(3);
  gfitR2->SetLineWidth(3);
  gfitR2->SetFillColor(3);
  gfitR2->SetFillStyle(3001);
  // gfitR2->Draw("p same");
  
  leg->AddEntry(gphotoR, "Photogrammetry phase-centres", "p");
  leg->AddEntry(gfitR, "Fitted phase-centres H-POL", "p");
  // leg->AddEntry(gfitR2, "Fitted phase-centres V-POL", "p");

  leg->Draw();

  c1->Print(Form("FittedValues_R_4steps.png"));
  c1->Print(Form("FittedValues_R_4steps.pdf"));


  TGraph *gphotoZ = new TGraph(48, antArray, z);
  gphotoZ->SetMarkerStyle(22);
  
  TGraphErrors *gfitZ = new TGraphErrors(48, antArray, zFit, antErr, deltaZErr);
  gfitZ->SetMarkerStyle(22);
  gfitZ->SetMarkerColor(kRed);

  TGraphErrors *gfitZ2 = new TGraphErrors(48, antArray, zFit2, antErr, deltaZErr2);
  gfitZ2->SetMarkerStyle(22);
  gfitZ2->SetMarkerColor(kBlue);


  gphotoZ->SetTitle(";Antenna;z [m]");
  
  gphotoZ->Draw("Ap");
  gphotoZ->SetMinimum(-4);
  gphotoZ->SetMaximum(4);
  gphotoZ->Draw("Ap");

  gfitZ->SetLineColor(2);
  gfitZ->SetLineWidth(2);
  gfitZ->SetFillColor(2);
  gfitZ->SetFillStyle(3001);
  gfitZ->Draw("p same");

  gfitZ2->SetLineColor(3);
  gfitZ2->SetLineWidth(3);
  gfitZ2->SetFillColor(3);
  gfitZ2->SetFillStyle(3001);
  // gfitZ2->Draw("p same");
  
  leg->Draw();

  c1->Print(Form("FittedValues_Z_4steps.png"));
  c1->Print(Form("FittedValues_Z_4steps.pdf"));


  TGraph *gphotoPHI = new TGraph(48, antArray, phi);
  gphotoPHI->SetMarkerStyle(22);
  
  TGraphErrors *gfitPHI = new TGraphErrors(48, antArray, phiFit, antErr, deltaPhiErr);
  gfitPHI->SetMarkerStyle(22);
  gfitPHI->SetMarkerColor(kRed);

  TGraphErrors *gfitPHI2 = new TGraphErrors(48, antArray, phiFit2, antErr, deltaPhiErr2);
  gfitPHI2->SetMarkerStyle(22);
  gfitPHI2->SetMarkerColor(kBlue);

  gphotoPHI->SetTitle(";Antenna;phi [rad]");
  
  gphotoPHI->Draw("Ap");
  gphotoPHI->SetMinimum(-2);
  gphotoPHI->SetMaximum(TMath::Pi()*2.1);
  gphotoPHI->Draw("Ap");

  gfitPHI->SetLineColor(2);
  gfitPHI->SetLineWidth(2);
  gfitPHI->SetFillColor(2);
  gfitPHI->SetFillStyle(3001);
  gfitPHI->Draw("p same");

  gfitPHI2->SetLineColor(3);
  gfitPHI2->SetLineWidth(3);
  gfitPHI2->SetFillColor(3);
  gfitPHI2->SetFillStyle(3001);
  // gfitPHI2->Draw("p same");
  
  leg->Draw();

  c1->Print(Form("FittedValues_PHI_4steps.png"));
  c1->Print(Form("FittedValues_PHI_4steps.pdf"));


  TGraph *gphotoCableDelay = new TGraph(48, antArray, antErr);
  gphotoCableDelay->SetMarkerStyle(22);
  
  TGraphErrors *gfitCableDelay = new TGraphErrors(48, antArray, deltaCableDelays , antErr, deltaCableDelaysErr);
  gfitCableDelay->SetMarkerStyle(22);
  gfitCableDelay->SetMarkerColor(kRed);

  TGraphErrors *gfitCableDelay2 = new TGraphErrors(48, antArray, deltaCableDelays2 , antErr, deltaCableDelaysErr2);
  gfitCableDelay2->SetMarkerStyle(22);
  gfitCableDelay2->SetMarkerColor(kBlue);

  gphotoCableDelay->SetTitle(";Antenna;Time offset [ns]");
  
  gphotoCableDelay->Draw("Ap");
  gphotoCableDelay->SetMinimum(-1);
  gphotoCableDelay->SetMaximum(1);
  gphotoCableDelay->Draw("Ap");

  gfitCableDelay->SetLineColor(2);
  gfitCableDelay->SetLineWidth(2);
  gfitCableDelay->SetFillColor(2);
  gfitCableDelay->SetFillStyle(3001);
  gfitCableDelay->Draw("p same");

  gfitCableDelay2->SetLineColor(3);
  gfitCableDelay2->SetLineWidth(3);
  gfitCableDelay2->SetFillColor(3);
  gfitCableDelay2->SetFillStyle(3001);
  // gfitCableDelay2->Draw("p same");
  
  leg->Draw();

  c1->Print(Form("FittedValues_CableDelay_4steps.png"));
  c1->Print(Form("FittedValues_CableDelay_4steps.pdf"));




  }

