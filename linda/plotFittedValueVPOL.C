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

void plotFittedValueVPOL();

void plotFittedValueVPOL(){


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

  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;

   for (unsigned int ant=0; ant<MAX_ANTENNAS; ++ant){
     antArray[ant] = (ant+1)*1.;
     // myGeomTool->getAntXYZ(ant, x[ant], y[ant], z[ant], pol);
     
     // r[ant] = myGeomTool->getAntR(ant, pol);
     // phi[ant] = myGeomTool->getAntPhiPosition(ant, pol);

   }

  TLegend *leg = new TLegend(0.6, 0.11, 0.89, 0.35);
  leg->SetFillColor(kWhite);

 
  TGraph *gphotoR = new TGraph(48, antArray, r);
  gphotoR->SetMarkerStyle(22);
  TGraph *gphotoZ = new TGraph(48, antArray, z);
  gphotoZ->SetMarkerStyle(22);
  TGraph *gphotoPHI = new TGraph(48, antArray, phi);
  gphotoPHI->SetMarkerStyle(22);
  TGraph *gphotoCableDelay = new TGraph(48, antArray, cableDelays);
  gphotoCableDelay->SetMarkerStyle(22);


  string filename[5]={"newLindaNumbers_4steps_VPOL_discone_2015_10_20_time_16_39_39.txt", "newLindaNumbers_4steps_VPOL_seavey_6kV_2015_10_20_time_14_29_34.txt", "newLindaNumbers_4steps_VPOL_seavey_10kV_2015_10_20_time_14_59_16.txt", "newLindaNumbers_4steps_VPOL_2015_10_15_time_10_51_50.txt", "newLindaNumbers_4steps_VPOL_ALL_2015_10_21_time_12_19_26.txt" };

  string which[5] = {"discone", "seavey6kV", "seavey10kV", "WAIS", "ALL"};



  
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

  for (int i=0;i<5;i++){
    
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
  c1->Print("FittedValues_VPOL_R_4steps.png");
  c1->Print("FittedValues_VPOL_R_4steps.pdf");

  c2->cd();
  leg->Draw();
  c2->Print("FittedValues_VPOL_Z_4steps.png");
  c2->Print("FittedValues_VPOL_Z_4steps.pdf");

  c3->cd();
  leg->Draw();
  c3->Print("FittedValues_VPOL_PHI_4steps.png");
  c3->Print("FittedValues_VPOL_PHI_4steps.pdf");

  c4->cd();
  leg->Draw();
  c4->Print("FittedValues_VPOL_CableDelay_4steps.png");
  c4->Print("FittedValues_VPOL_CableDelay_4steps.pdf");


  }

