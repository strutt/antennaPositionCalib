void plotAngularResolution(){


   TChain *resChain = new TChain("angResTree");

   char resName[500];
   char fName[200];

   //   int sec1[17] = {04, 04, 04, 04, 04, 04, 10, 10, 12, 12, 12, 14, 14, 14, 16, 16, 17};
   //   int sec2[17] = {51, 51, 57, 57, 57, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59};
   int temp = -1;
   //   int sec1[25] = {24, 30, 30, 25, 25, 25, 26, 26, 25, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 29, 30, 30, 31};

   //   int sec1[17]={36, 36, 36, 36, 36, 33, 34, 34, 35, 36, 35, 36, 36, 37, 37, 37, 37};
   int sec1[3]= {19,21,21};

   for (int run=331;run<355;++run){ // 331-355 145-161 151-154

     temp = run - 151;
     //     cout << sec1[temp] << endl;
     //     sprintf(fName,"VPOL_10kVSeavey_2015_11_05_time_19_47_53");//  VPOL_10kVSeavey_2015_11_06_time_10_42_48
     // sprintf(fName,"VPOL_10kVSeavey_2015_11_06_time_10_42_48");
     // sprintf(fName,"LDBHPOL_10kVSeavey_2015_11_19_time_15_30_45");
     //     sprintf(fName,"VPOL_10kVSeavey_2015_11_19_time_11_49_04");
     sprintf(fName,"HPOL_seavey_pulses_2015_11_25");
     
          sprintf(resName,"/unix/anita3/treesForLinda/generateAngularResolutionTree_run%d-%dPlots.root",run,run);
     //     sprintf(resName,"/unix/anita3/strutt/antennaPositionCalibPlots/ldbPulses/%s/generateAngularResolutionTreeVPOLPlots_%d_2015-11-09_14-36-%02d.root",fName,run, sec2[temp]);
     //     sprintf(resName,"/unix/anita3/strutt/antennaPositionCalibPlots/newFittingOrder/newLindaNumbers_4steps_%s/generateAngularResolutionTreeVPOLPlots_%d_2015-11-23_12-42-%02d.root",fName,run, sec1[temp]);
     //sprintf(resName,"/unix/anita3/strutt/antennaPositionCalibPlots/ldbPulses/%s/generateAngularResolutionTreeHPolLDBPlots_%d_2015-11-25_12-25-%02d.root",fName,run, sec1[temp]);

     cout << resName << endl;
     resChain->Add(resName);
     
  }


   TH1D *hPhi = new TH1D("hPhi", "", 400, -2, 2);
   TH1D *hTheta = new TH1D("hTheta", "", 400, -2, 2);

   resChain->Draw("deltaPhiDeg >> hPhi", "", "off");
   resChain->Draw("deltaThetaDeg >> hTheta", "", "off");
   

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);
   TCanvas *c1 = new TCanvas("c1");
   hPhi->SetTitle("#phi_{expected}-#phi_{zoom};#delta #phi (degrees);Events per bin");
   hPhi->Draw("");
   hPhi->Fit("gaus", "", "", -1, 1);
   c1->Print(Form("AngularResolution_phi_4steps_%s.png", fName));
   c1->Print(Form("AngularResolution_phi_4steps_%s.pdf", fName));


   hTheta->SetTitle("#theta_{expected}-#theta_{zoom};#delta #theta (degrees);Events per bin");
   hTheta->Draw("");
   hTheta->Fit("gaus");
   c1->Print(Form("AngularResolution_theta_4steps_%s.png", fName));
   c1->Print(Form("AngularResolution_theta_4steps_%s.pdf", fName));



}
