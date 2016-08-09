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

void plotFittedValueHPOL2();

void plotFittedValueHPOL2(){


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
     // phi[ant] = myGeomTool->getAntPhiPositionRelToAftFore(ant, pol);

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

  const int numfile =3;
  //  string filename[numfile]={"newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_05_time_19_47_53.txt", "newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_06_time_10_42_48.txt", "newLindaNumbers_4steps_zDisplaced_VPOL_10kVSeavey_2015_11_11_time_15_38_30.txt", "newLindaNumbers_4steps_2015_10_13_time_14_30_54.txt"};
  // string filename[numfile]={"newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_19_time_11_49_04.txt", "final/newLindaNumbers_4steps_HPOL_2015_11_13_time_17_19_58.txt", "newLindaNumbers_4steps_HPOL_2015_11_19_time_15_06_17.txt", "newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_10_16_23.txt", "newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45.txt"};
  string filename[numfile]={"final/newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_19_time_11_49_04.txt", "final/newLindaNumbers_4steps_HPOL_2015_11_19_time_15_06_17.txt", "final/newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45.txt"};
    
  //  string filename[numfile] = {"final/newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_19_time_11_49_04.txt", "newFittingMethod/newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_27_time_11_49_56.txt", "newLindaNumbers_4steps_HPOL_2015_11_19_time_15_06_17.txt", "newLindaNumbers_4steps_WAISHPOL_2015_11_27_time_12_13_09.txt", "newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45.txt", "newFittingMethod/newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_27_time_11_50_07.txt", "experiment/newLindaNumbers_LDBHPOL_2015_11_27_time_17_03_06.txt"};

  string which[numfile] = {"FIT V-POL", "FIT WAIS H-POL", "FIT LDB H-POL"};

  double diff[48];


  TCanvas *c1 = new TCanvas("c1", "", 1200, 750);
  gphotoR->SetTitle(";Antenna;#Delta r [m]");  
  gphotoR->Draw("Ap");
  gphotoR->SetMinimum(-0.03);
  gphotoR->SetMaximum(+0.03);
  gphotoR->Draw("Ap");

  TCanvas *c2 = new TCanvas("c2", "", 1200, 750);
  gphotoZ->SetTitle(";Antenna;#Delta z [m]");
  gphotoZ->Draw("Ap");
  gphotoZ->SetMinimum(-0.01);
  gphotoZ->SetMaximum(+0.01);
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

  Int_t colors[10];//{kRed, kOrange, kGreen, kBlue, kCyan, kGray};
  for (int i=0;i<10;i++) colors[i] = 52+i*8;
    

  for (int ifile=0;ifile<numfile;ifile++){
    
    ifstream myInput(Form("/home/lindac/ANITA/Software/EventCorrelator/macros/lindaMacros/%s", filename[ifile].c_str()));
    
    if (myInput.is_open()){
      for (int i=0;i<48;i++){
	myInput >> antArray[i] >> deltaR[i] >> deltaZ[i] >> deltaPhi[i] >> deltaCableDelays[i];
	if (ifile==2) cout << antArray[i] << " " << deltaR[i]+r[i] << "\t" << deltaZ[i]+z[i] << "\t" << (deltaPhi[i]+phi[i])*TMath::RadToDeg() << endl;

	if (ifile==1) diff[i] = deltaPhi[i];
	else if (ifile==2) diff[i]-=deltaPhi[i];
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
   gfitR->SetMarkerColor(colors[ifile]);

   TGraph *gfitZ = new TGraph(48, antArray, zFit);
   gfitZ->SetMarkerStyle(22);
   gfitZ->SetMarkerColor(colors[ifile]);

   TGraph *gfitPHI = new TGraph(48, antArray, phiFit);
   gfitPHI->SetMarkerStyle(22);
   gfitPHI->SetMarkerColor(colors[ifile]);

   TGraph *gfitT = new TGraph(48, antArray, deltaCableDelays);
   gfitT->SetMarkerStyle(22);
   gfitT->SetMarkerColor(colors[ifile]);

   leg->AddEntry(gfitR, Form("Fitted w/ %s", which[ifile].c_str()), "p");

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
  c1->Print("FittedValues_VPOL_R_4steps_compare.png");
  c1->Print("FittedValues_VPOL_R_4steps_compare.pdf");

  c2->cd();
  leg->Draw();
  c2->Print("FittedValues_VPOL_Z_4steps_compare.png");
  c2->Print("FittedValues_VPOL_Z_4steps_compare.pdf");

  c3->cd();
  leg->Draw();
  c3->Print("FittedValues_VPOL_PHI_4steps_compare.png");
  c3->Print("FittedValues_VPOL_PHI_4steps_compare.pdf");

  c4->cd();
  leg->Draw();
  c4->Print("FittedValues_VPOL_CableDelay_4steps_compare.png");
  c4->Print("FittedValues_VPOL_CableDelay_4steps_compare.pdf");


  TCanvas *c5 = new TCanvas("c5");

   TGraph *gfitPHIDIFF = new TGraph(48, antArray, diff);
   gfitPHIDIFF->SetMarkerStyle(22);
   gfitPHIDIFF->SetMarkerColor(kBlue);

  gfitPHIDIFF->SetMinimum(-0.1);
  gfitPHIDIFF->SetMaximum(+0.1);

  double mean = 0;
  for (int i=0;i<48;i++) mean+=diff[i]/48.;

  cout << "Mean diff : " << mean*TMath::RadToDeg() << " degree" << endl;

  gfitPHIDIFF->Draw("Ap");
  c5->Print("diff.png");


  }

