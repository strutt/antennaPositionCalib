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

void plotFittedValueHPOL();

void plotFittedValueHPOL(){


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

  Double_t deltaR[MAX_ANTENNAS]={0};  
  Double_t deltaZ[MAX_ANTENNAS]={0};  
  Double_t deltaPhi[MAX_ANTENNAS]={0};  
  Double_t deltaCableDelays[MAX_ANTENNAS]={0};

  Double_t antArray[MAX_ANTENNAS] = {0.};
  Double_t cableDelays[MAX_ANTENNAS]={0};

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

   for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
     antArray[ant] = (ant+1)*1.;
     // myGeomTool->getAntXYZ(ant, x[ant], y[ant], z[ant], pol);
     
     // r[ant] = myGeomTool->getAntR(ant, pol);
     // phi[ant] = myGeomTool->getAntPhiPosition(ant, pol);

   }

  TLegend *leg = new TLegend(0.6, 0.75, 0.89, 0.89);
  leg->SetFillColor(kWhite);

 
  TGraph *gphotoR = new TGraph(48, antArray, r);
  gphotoR->SetMarkerStyle(22);
  TGraph *gphotoZ = new TGraph(48, antArray, z);
  gphotoZ->SetMarkerStyle(22);
  TGraph *gphotoPHI = new TGraph(48, antArray, phi);
  gphotoPHI->SetMarkerStyle(22);
  TGraph *gphotoCableDelay = new TGraph(48, antArray, cableDelays);
  gphotoCableDelay->SetMarkerStyle(22);

  const int numfile = 2;

  //  string filename[numfile]={"newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_05_time_19_47_53.txt", "newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_06_time_10_42_48.txt", "newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_11_time_15_38_30.txt", "newLindaNumbers_4steps_2015_10_13_time_14_30_54.txt"};
  string filename[numfile]={"final/newLindaNumbers_4steps_HPOL_2015_11_13_time_17_19_58.txt", "newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_10_16_23.txt"};

  string which[numfile] = {"FIT WAIS H-POL", "FIT LDB H-POL"};

  
  TCanvas *c1 = new TCanvas("c1", "", 1200, 750);
  gphotoR->SetTitle(";Antenna;#Delta r [m]");  
  gphotoR->Draw("Ap");
  gphotoR->SetMinimum(-0.2);
  gphotoR->SetMaximum(+0.2);
  gphotoR->Draw("Ap");

  TCanvas *c2 = new TCanvas("c2", "", 1200, 750);
  gphotoZ->SetTitle(";Antenna;#Delta z [m]");
  gphotoZ->Draw("Ap");
  gphotoZ->SetMinimum(-0.2);
  gphotoZ->SetMaximum(+0.2);
  gphotoZ->Draw("Ap");

  TCanvas *c3 = new TCanvas("c3", "", 1200, 750);
  gphotoPHI->SetTitle(";Antenna;#Delta phi [rad]");  
  gphotoPHI->Draw("Ap");
  gphotoPHI->SetMinimum(-0.2);
  gphotoPHI->SetMaximum(+0.2);
  gphotoPHI->Draw("Ap");

  TCanvas *c4 = new TCanvas("c4", "", 1200, 750);
  gphotoCableDelay->SetTitle(";Antenna;Time offset [ns]");
  gphotoCableDelay->Draw("Ap");
  gphotoCableDelay->SetMinimum(-0.5);
  gphotoCableDelay->SetMaximum(+0.5);
  gphotoCableDelay->Draw("Ap");

  leg->AddEntry(gphotoR, "Photogrammetry", "p");

  Color_t colors[5]={kRed, kOrange, kGreen, kBlue, kCyan};

  for (int i=0;i<numfile;i++){
    
    ifstream myInput(Form("/home/lindac/ANITA/Software/EventCorrelator/macros/lindaMacros/%s", filename[i].c_str()));
    
    if (myInput.is_open()){
      for (int i=0;i<48;i++){
	myInput >> antArray[i] >> deltaR[i] >> deltaZ[i] >> deltaPhi[i] >> deltaCableDelays[i];
	cout << antArray[i] << " " << deltaR[i] << endl;
      }
    }
    
   for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
     
    rFit[ant] = r[ant] + deltaR[ant];
    zFit[ant] = z[ant] + deltaZ[ant];
    phiFit[ant] = phi[ant] + deltaPhi[ant];
    antArray[ant] = (ant+1)*1.;
    
    
   }
    
   TGraph *gfitR = new TGraph(48, antArray, rFit);
   gfitR->SetMarkerStyle(22);
   gfitR->SetMarkerColor(colors[i]);

   TGraph *gfitZ = new TGraph(48, antArray, zFit);
   gfitZ->SetMarkerStyle(22);
   gfitZ->SetMarkerColor(colors[i]);

   TGraph *gfitPHI = new TGraph(48, antArray, phiFit);
   gfitPHI->SetMarkerStyle(22);
   gfitPHI->SetMarkerColor(colors[i]);

   TGraph *gfitT = new TGraph(48, antArray, deltaCableDelays);
   gfitT->SetMarkerStyle(22);
   gfitT->SetMarkerColor(colors[i]);

   leg->AddEntry(gfitR, Form("Fitted w/ %s", which[i].c_str()), "p");

   c1->cd();
   gfitR->Draw("p same");

   c2->cd();
   gfitZ->Draw("p same");

   c3->cd();
   gfitPHI->Draw("p same");

   c4->cd();
   gfitT->Draw("p same");

   
  }

  c1->cd();
  leg->Draw();
  c1->Print("FittedValues_HPOL_R_4steps_compare.png");
  c1->Print("FittedValues_HPOL_R_4steps_compare.pdf");

  c2->cd();
  leg->Draw();
  c2->Print("FittedValues_HPOL_Z_4steps_compare.png");
  c2->Print("FittedValues_HPOL_Z_4steps_compare.pdf");

  c3->cd();
  leg->Draw();
  c3->Print("FittedValues_HPOL_PHI_4steps_compare.png");
  c3->Print("FittedValues_HPOL_PHI_4steps_compare.pdf");

  c4->cd();
  leg->Draw();
  c4->Print("FittedValues_HPOL_CableDelay_4steps_compare.png");
  c4->Print("FittedValues_HPOL_CableDelay_4steps_compare.pdf");


  }

