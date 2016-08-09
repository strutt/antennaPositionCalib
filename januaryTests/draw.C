void draw(){

  // const int numDirs = 5;
  // TString dirNames[numDirs] = {"photogrammetryNumbers", "frameworkNumbers", "newLindaNumbers_4steps_WAISHPOL_2015_12_17_time_18_22_28", "newLindaNumbers_4steps_VPOL_10kVSeavey_2015_12_17_time_17_41_28", "newLindaNumbers_LDBHPOL_2015_12_17_time_17_38_21"};
  // TString titles[numDirs] = {"Photogrammetry", "Framework", "New WAIS", "New 10kV VPOL", "New LDB HPol"};  
  // const int numDirs = 2;
  // TString dirNames[numDirs] = {"frameworkNumbers", "newCoordinatesLDB/frameworkNumbers"};
  // TString titles[numDirs] = {"Framework_old_LDB_location", "Framework_new_LDB_location"};

  // const int numDirs = 10;  
  
  // TString dirNames[numDirs]={"frameworkNumbers",
  // 			     "newLindaNumbers_4steps_2015_10_13_time_14_30_54",
  // 			     "newLindaNumbers_LDBHPOL_2015_12_17_time_17_38_21",
  // 			     "newLindaNumbers_4steps_WAISHPOL_2015_12_17_time_18_22_28",
  // 			     "newLindaNumbers_4steps_WAISHPOL_2016_01_11_time_16_01_54",
  // 			     "newLindaNumbers_4steps_WAISHPOL_NEW3_2016_01_11_time_19_30_31",
  // 			     "newLindaNumbers_4steps_WAISHPOL_NEW8_2016_01_18_time_17_10_01",
  // 			     "newLindaNumbers_4steps_WAISHPOL_NEW7_2016_01_18_time_18_18_35",
  // 			     "newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11",
  // 			     "cosminGpsTrees/newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11"
  // };

  // const int numDirs = 2;  

  // TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11",
  // 			     "cosminGpsTrees/newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11"};  
  
  // 			     // "newLindaNumbers_4steps_2015_10_13_time_14_30_54",
  // 			     // "newLindaNumbers_4steps_VPOL_10kVSeavey_TEEEEEEEST_2016_01_11_time_14_56_23",
  // 			     // "newLindaNumbers_4steps_WAISHPOL_2015_12_17_time_18_22_28",
  // 			     // "newLindaNumbers_LDBHPOL_2016_01_11_time_14_42_08"};
  // TString titles[numDirs] = {"Framework", "FrameworkLindasDts","HPolLDB", "Wais", "JanuaryWAIS", "WAIS_NEW3", "WAIS_NEW8", "WAIS_NEW7", "WAIS_NEW10", "cosminTreesWAIS_NEW10"};//, "VPolLDB", "Wais", "HPolLDB"};


  // // TString titles[numDirs] = {"WAIS_NEW10", "cosminTreesWAIS_NEW10"};//, "VPolLDB", "Wais", "HPolLDB"};

  // const int numDirs = 3;
  // // TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11",
  // // 			     "newLindaNumbers_LDBHPOL_NEW10_2016_01_19_time_20_38_29",
  // // 			     "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_2016_01_19_time_20_39_33"};
  // TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV2_2016_01_21_time_14_53_41",
  // 			     "newLindaNumbers_LDBHPOL_NEW10_cosminV2_2016_01_21_time_14_03_25",
  // 			     "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_cosminV2_2016_01_21_time_14_30_18"};

  // TString titles[numDirs] = {"WAIS_NEW10_cosminV2", "LDBHPOL_NEW10_cosminV2", "10kVSeavey_NEW10_cosminV2"};


  // const int numDirs = 3;
  // // TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11",
  // // 			     "newLindaNumbers_LDBHPOL_NEW10_2016_01_19_time_20_38_29",
  // // 			     "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_2016_01_19_time_20_39_33"};
  // TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV3_2016_01_25_time_12_24_25",
  // 			     "newLindaNumbers_LDBHPOL_NEW10_cosminV3_2016_01_25_time_11_33_16",
  // 			     "newLindaNumbers_4steps_VPOL_10kVSeavey_NEW10_cosminV3_2016_01_25_time_11_44_07"};

  // TString titles[numDirs] = {"WAIS_NEW10_cosminV3", "LDBHPOL_NEW10_cosminV3", "10kVSeavey_NEW10_cosminV3"};

  const int numDirs = 1;
  TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15"};

  TString titles[numDirs] = {"Bug fixes"};
  
  for(int dirInd=0; dirInd < numDirs; dirInd++){
    auto cWais = new TChain("angResTree");
    auto cHPolLDB = new TChain("angResTree");
    auto cVPol = new TChain("angResTree");

    // cWais->Add("frameworkNumbers/generateAngularResolutionTreePlots_*");
    // cHPolLDB->Add("frameworkNumbers/generateAngularResolutionTreeHPolLDBPlots_*");  
    // cVPol->Add("frameworkNumbers/generateAngularResolutionTreeVPOLPlots_*");

    // cWais->Add("photogrammetryNumbers/generateAngularResolutionTreePlots_*");
    // cHPolLDB->Add("photogrammetryNumbers/generateAngularResolutionTreeHPolLDBPlots_*");    
    // cVPol->Add("photogrammetryNumbers/generateAngularResolutionTreeVPOLPlots_*");

    // cWais->Add(dirNames[dirInd] + "/generateAngularResolutionTreePlots_*");
    // cHPolLDB->Add(dirNames[dirInd] + "/generateAngularResolutionTreeHPolLDBPlots_*");
    // cVPol->Add(dirNames[dirInd] + "/generateAngularResolutionTreeVPOLPlots_*");

    // cWais->Add("newCoordinatesLDB/"+dirNames[dirInd] + "/generateAngularResolutionTreePlots_*");
    // cHPolLDB->Add("newCoordinatesLDB/"+dirNames[dirInd] + "/generateAngularResolutionTreeHPolLDBPlots_*");
    // cVPol->Add("newCoordinatesLDB/"+dirNames[dirInd] + "/generateAngularResolutionTreeVPOLPlots_*");
    cWais->Add(dirNames[dirInd] + "/generateAngularResolutionTreePlots_*");
    cHPolLDB->Add(dirNames[dirInd] + "/generateAngularResolutionTreeHPolLDBPlots_*");
    cVPol->Add(dirNames[dirInd] + "/generateAngularResolutionTreeVPOLPlots_*");

    cout << dirNames[dirInd].Data() << "\t" << cWais->GetEntries() << "\t" << cHPolLDB->GetEntries() << "\t" << cVPol->GetEntries() << endl;

    const double angleRange = 5; //3; //180;
    const int numBinsY = 128*1.6666; //128*60; //64;
    // const int numBinsX = 128; //64;
    const int numBinsX = 128; //64;        

    TH2D* hWais = new TH2D("hWais"+titles[dirInd], "WAIS pulses", numBinsX, 0, cWais->GetEntries(), numBinsY, -angleRange, angleRange);
    hWais->GetXaxis()->SetTitle("Entry");
    hWais->GetYaxis()->SetTitle("#delta#phi (Degrees)");    
    // TH2D* hWais = new TH2D("hWais"+titles[dirInd], "WAIS pulses", numBinsX, 0, 360, numBinsY, -angleRange, angleRange);    
    // cWais->Draw("deltaPhiDeg>>hWais"+titles[dirInd], "TMath::Abs(deltaPhiDeg) < 3 && heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153 && mrms < 0.1 && mrms > 0 && brms < 40e-3", "goff");
    // cWais->Draw("deltaPhiDeg:zoomPhiDeg>>hWais"+titles[dirInd], "heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153", "goff");
    // cWais->Draw("deltaPhiDeg:Entry$>>hWais"+titles[dirInd], "heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153", "goff");
    cWais->Draw("triggeredPhiDeg-zoomPhiDeg:Entry$>>hWais"+titles[dirInd], "heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153", "goff");     

    TH2D* hHPolLDB = new TH2D("hHPolLDB"+titles[dirInd], "HPol LDB", numBinsX, 0, cHPolLDB->GetEntries(), numBinsY, -angleRange, angleRange);
    // TH2D* hHPolLDB = new TH2D("hHPolLDB"+titles[dirInd], "HPol LDB", numBinsX, 0, 360, numBinsY, -angleRange, angleRange);    
    hHPolLDB->GetXaxis()->SetTitle("Entry");
    hHPolLDB->GetYaxis()->SetTitle("#delta#phi (Degrees)");    

    // cHPolLDB->Draw("deltaPhiDeg>>hHPolLDB"+titles[dirInd], "TMath::Abs(deltaPhiDeg) < 3 && heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153 && mrms < 0.1 && mrms > 0 && brms < 40e-3", "goff");
    cHPolLDB->Draw("deltaPhiDeg:Entry$>>hHPolLDB"+titles[dirInd], "heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153", "goff");
    // cHPolLDB->Draw("deltaPhiDeg:zoomPhiDeg>>hHPolLDB"+titles[dirInd], "heading > -50 && l3TrigPatternH > 0 && run!=150 && run!=153", "goff");    
    hHPolLDB->SetLineColor(kRed);

    TH2D* hVPol = new TH2D("hVPol"+titles[dirInd], "VPol LDB", numBinsX, 0, cVPol->GetEntries(), numBinsY, -angleRange, angleRange);
    // TH2D* hVPol = new TH2D("hVPol"+titles[dirInd], "VPol LDB", numBinsX, 0, 360, numBinsY, -angleRange, angleRange);
    hVPol->GetXaxis()->SetTitle("Entry");
    hVPol->GetYaxis()->SetTitle("#delta#phi (Degrees)");        
    // cVPol->Draw("deltaPhiDeg>>hVPol"+titles[dirInd], "TMath::Abs(deltaPhiDeg) < 3 && heading > -50 && l3TrigPattern > 0 && run!=150 && run!=153 && mrms < 0.1 && mrms > 0 && brms < 40e-3", "goff");
    // cVPol->Draw("deltaPhiDeg:zoomPhiDeg>>hVPol"+titles[dirInd], "heading > -50 && l3TrigPattern > 0 && run!=150 && run!=153", "goff");
    cVPol->Draw("deltaPhiDeg:Entry$>>hVPol"+titles[dirInd], "heading > -50 && l3TrigPattern > 0 && run!=150 && run!=153", "goff");    
    hVPol->SetLineColor(kBlue);
    

    const int nHists = 3;
    TH1D* hFramework[nHists] = {hHPolLDB->ProjectionY(), hWais->ProjectionY(), hVPol->ProjectionY()};
    for(int hInd=0; hInd < nHists; hInd++){
      if(hFramework[hInd]->Integral() > 0){
	hFramework[hInd]->Scale(1./hFramework[hInd]->Integral());
      }
    }
    
    // const int nHists = 2;
    // TH1D* hFramework[nHists] = {hHPolLDB, hWais};
    
    auto c1 = RootTools::drawHistsWithStatsBoxes(nHists, hFramework, "", "mre");
    c1->Update();
    TLegend* l = c1->BuildLegend();

    // hFramework[0]->SetTitle(titles[dirInd]+" Numbers; #delta#phi (Degrees); Fraction of events / bin");
    hFramework[0]->SetTitle(titles[dirInd]+" Numbers; #delta#phi (Degrees); Fraction of events / bin");

    
    l->Draw();
    c1->Update();
  }  
}
