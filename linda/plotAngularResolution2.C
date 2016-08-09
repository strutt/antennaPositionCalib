//////////////////////////////////////////////
// arSample = angular resolution sample 
// pcSample = phase centres sample
/////////////////////////////////////////////

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TPaveStats.h"
#include "TChain.h"
#include "TList.h"
#include "TF1.h"

#include <iostream>
using namespace std;

void plotAngularResolution_histos(string arSample, string pcSample , TH1D *hPhi, TH1D *hTheta);

void plotAngularResolution_compare(string AngResSample, string angle);
void plotAngularResolution_VPOL(string AngResSample, string angle);

void plotAngularResolution2(){

  plotAngularResolution_compare("WAIS","phi");
  plotAngularResolution_compare("WAIS","theta");
  plotAngularResolution_compare("LDBHPOL","phi");
  plotAngularResolution_compare("LDBHPOL","theta");
  plotAngularResolution_VPOL("LDBVPOL", "phi");
  plotAngularResolution_VPOL("LDBVPOL", "theta");
}

void plotAngularResolution_VPOL(string AngResSample, string angle){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  
  // Looking at WAIS pulses with WAIS and LDB HPOL phase-centres

  TH1D *hPhiLDB = new TH1D("hPhiLDB", "", 100, -2, 2);
  TH1D *hThetaLDB = new TH1D("hThetaLDB", "", 200, -2, 2);

  plotAngularResolution_histos(AngResSample, "LDBVPOL", hPhiLDB, hThetaLDB);
  hPhiLDB->SetName("hPhiLDB");
  hThetaLDB->SetName("hThetaLDB");
  

  if (angle=="theta"){
    TH1D* hLDB = (TH1D*)hThetaLDB->Clone();
  } else {
    TH1D* hLDB = (TH1D*)hPhiLDB->Clone();

  }

  hLDB->SetName("hLDB");

 
  TF1 *f1 = new TF1("f1", "gaus", -2, 2);
  TF1 *f2 = new TF1("f2", "gaus", -2, 2);

  TCanvas *c1 = new TCanvas("c1");
  hLDB->SetTitle(Form("#%s_{expected}-#%s_{zoom}: %s pulses;#delta #phi (degrees);Events per bin", angle.c_str(), angle.c_str(), AngResSample.c_str()));
  hLDB->SetLineColor(kBlue);
  hLDB->SetLineWidth(2);
  hLDB->Draw("");
  f1->SetLineColor(kBlue);
  hLDB->Fit(f1, "", "", -1, 1);

  gPad->Update();

  TPaveText *pt1 = new TPaveText(.13,.6,.4,.89, "NDC");
  pt1->SetTextColor(kBlue);
  writept(f1, pt1, "LDB");
  pt1->Draw("same");


    c1->Print(Form("angularResolution/AngularResolution_%ssample_%s.png", AngResSample.c_str(), angle.c_str()));
    c1->Print(Form("angularResolution/AngularResolution_%ssample_%s.pdf", AngResSample.c_str(), angle.c_str()));

}


void plotAngularResolution_compare(string AngResSample, string angle){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  
  // Looking at WAIS pulses with WAIS and LDB HPOL phase-centres

  TH1D *hPhiWAIS = new TH1D("hPhiWAIS", "", 100, -2, 2);
  TH1D *hThetaWAIS = new TH1D("hThetaWAIS", "", 200, -2, 2);
  TH1D *hPhiLDB = new TH1D("hPhiLDB", "", 100, -2, 2);
  TH1D *hThetaLDB = new TH1D("hThetaLDB", "", 200, -2, 2);

  plotAngularResolution_histos(AngResSample, "WAIS",  hPhiWAIS, hThetaWAIS);
  hPhiWAIS->SetName("hPhiWAIS");
  hThetaWAIS->SetName("hThetaWAIS");

  plotAngularResolution_histos(AngResSample, "LDBHPOL", hPhiLDB, hThetaLDB);
  hPhiLDB->SetName("hPhiLDB");
  hThetaLDB->SetName("hThetaLDB");
  

  if (angle=="theta"){
    TH1D* hWAIS = (TH1D*)hThetaWAIS->Clone();
    TH1D* hLDB = (TH1D*)hThetaLDB->Clone();
  } else {
    TH1D* hWAIS = (TH1D*)hPhiWAIS->Clone();
    TH1D* hLDB = (TH1D*)hPhiLDB->Clone();

  }

  hWAIS->SetName("hWAIS");
  hLDB->SetName("hLDB");

 
  TF1 *f1 = new TF1("f1", "gaus", -2, 2);
  TF1 *f2 = new TF1("f2", "gaus", -2, 2);

  TCanvas *c1 = new TCanvas("c1");
  hWAIS->SetTitle(Form("#%s_{expected}-#%s_{zoom}: %s pulses;#delta #phi (degrees);Events per bin", angle.c_str(), angle.c_str(), AngResSample.c_str()));
  hWAIS->SetLineColor(kBlue);
  hWAIS->SetLineWidth(2);
  hWAIS->Draw("");
  f1->SetLineColor(kBlue);
  hWAIS->Fit(f1, "", "", -1, 1);

  gPad->Update();

  TPaveText *pt1 = new TPaveText(.13,.6,.4,.89, "NDC");
  pt1->SetTextColor(kBlue);
  writept(f1, pt1, "WAIS");
  pt1->Draw("same");

  hLDB->SetTitle(Form("#%s_{expected}-#%s_{zoom}: %s pulses;#delta #phi (degrees);Events per bin", angle.c_str(), angle.c_str(), AngResSample.c_str()));
  hLDB->SetLineColor(kRed);
  hLDB->SetLineWidth(2);
  hLDB->Draw("same");
  hLDB->Fit(f2, "", "", -1, 1);
  f2->SetLineColor(kRed);
  gPad->Update();

  TPaveText *pt2 = new TPaveText(.13,.29,.4,.59, "NDC");
  pt2->SetTextColor(kRed);
  writept(f2, pt2, "LDB");
  pt2->Draw("same");
    

    c1->Print(Form("angularResolution/AngularResolution_%ssample_%s.png", AngResSample.c_str(), angle.c_str()));
    c1->Print(Form("angularResolution/AngularResolution_%ssample_%s.pdf", AngResSample.c_str(), angle.c_str()));

}

void plotAngularResolution_histos(string arSample, string pcSample, TH1D *hPhi, TH1D *hTheta){
  
  int firstRun, lastRun;
  string folderName;
  //  string baseName = "/unix/anita3/";
  string baseName = "/unix/anita3/strutt/antennaPositionCalibPlots/agreedNumbers/";

  // choose the phase-centres to use
  if (pcSample=="WAIS"){
    //    folderName = "strutt/antennaPositionCalibPlots//newFittingOrder/newLindaNumbers_4steps_WAISHPOL_2015_11_19_time_15_06_17";
    //    folderName = "newLindaNumbers_4steps_WAISHPOL_2015_11_27_time_12_13_09";
    //    folderName="strutt/antennaPositionCalibPlots/newFittingOrder/newLindaNumbers_4steps_HPOL_2015_11_13_time_17_19_58/";
    //    folderName="treesForLinda/";
    folderName="newLindaNumbers_4steps_WAISHPOL_2015_12_03_time_00_43_25";
  } else if (pcSample=="LDBHPOL") {
    //   folderName = "strutt/antennaPositionCalibPlots/newFittingOrder/newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45";
    //    folderName = "strutt/antennaPositionCalibPlots/ldbPulses/newLindaNumbers_4steps_LDBHPOL_10kVSeavey_2015_11_27_time_11_50_07";
    folderName="newLindaNumbers_LDBHPOL_2015_12_03_time_00_02_44";
  } else if (pcSample=="LDBVPOL") {
    //    folderName = "strutt/antennaPositionCalibPlots/newFittingOrder/newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_19_time_11_49_04";
    //  folderName = "antennaPositionCalibPlots/ldbPulses/newLindaNumbers_4steps_VPOL_10kVSeavey_2015_11_27_time_11_49_56/";
    folderName="newLindaNumbers_4steps_VPOL_10kVSeavey_2015_12_03_time_00_19_55";
  }

  hPhi->SetName("hPhi");
  hTheta->SetName("hTheta");


  // choose the sample for the angular resolution plots
  if (arSample=="WAIS"){
    firstRun = 331;
    lastRun = 354;
  } else if (arSample=="LDBHPOL") {
    firstRun = 151;
    lastRun = 153;
  } else if (arSample=="LDBVPOL") {
    firstRun = 145;
    lastRun = 161;
  }

  
  TChain *resChain = new TChain("angResTree");

  char resName[300];
  for (int irun=firstRun; irun<lastRun; irun++){
    sprintf(resName,"%s%s/*%d*root", baseName.c_str(), folderName.c_str(), irun);
    cout << resName << endl;
    resChain->Add(resName);
  }


  resChain->Draw("deltaPhiDeg >> hPhi", "", "off");
  resChain->Draw("deltaThetaDeg >> hTheta", "", "off");

  // TCanvas *c1 = new TCanvas ("c1");
  // resChain->Draw("deltaPhiDeg:eventNumber", "", "colz");
  // c1->Print("temp.png");


  resChain->Delete();

}


void writept(TF1 *f1, TPaveText *pt, string pc){

  pt->SetFillColor(kWhite);
  pt->AddText(Form("%s PHASE-CENTRES", pc.c_str()));
  pt->AddText(Form("#chi^2 / ndf %20.1f / %4d", f1->GetChisquare(), f1->GetNDF()));
  pt->AddText(Form("Constant %18.1f #pm %5.1f", f1->GetParameter(0), f1->GetParError(0)) );
  pt->AddText(Form("Mean %19.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)) );
  pt->AddText(Form("Width %18.4f #pm %5.4f", f1->GetParameter(2), f1->GetParError(2)) );

}
