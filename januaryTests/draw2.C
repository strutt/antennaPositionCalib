void draw2(){

  gStyle->SetOptStat("mre");
  
  const int numDirs = 4;
  TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15",
			     "testing/corInterp/latestNums",
			     "testing/corInterp/photogrammetry",
			     "../cosminsTrees"};

  // TString titles[numDirs] = {"bf1", "bins0.5", "bigBins", "corInterp", "corInterp_photo"};
  // TString titles[numDirs] = {"bf1", "corInterpLatest", "corInterp_photo", "cosmin"};
  TString titles[numDirs] = {"#deltat binning fix", "Correlation interpolation", "Correlation interpolation (photogrammetry)", "Cosmin's Reco"};


  std::vector<TH1D*> hists;

  // auto c1 = new TCanvas();
  for(int dirInd=0; dirInd < numDirs; dirInd++){
    auto cWais = new TChain("angResTree");
    cWais->Add(dirNames[dirInd] + "/generateAngularResolutionTreePlots_352*root");

    TString varName = "deltaPhiDeg";
    // TString varName = "triggeredPhiDeg - phiExpected";
    // TString varName = "phiExpected - triggeredPhiDeg";
    // TString varName = "phiExpected - zoomPhiDeg";
    TString cuts = "heading > -50 && l3TrigPatternH > 0";
    
    if(cWais->GetEntries() == 0){
      cWais = new TChain("pointed");
      cWais->Add(dirNames[dirInd] + "/352.root");      
      varName = "-WAIS.dphi";
      cuts = "heading > -50";
    }
    
    const double angleRange = 5; //3; //180;
    const int numBins = 50;

    // // new TCanvas();
    // TString name = TString::Format("h1_%d", dirInd);
    // TH1D* h1 = new TH1D(name,
    // 			"triggeredPhiDeg-zoomPhiDeg:triggeredPhiDeg",
    // 			// 128, 0, 360,
    // 			51, -5, 5);
    // cWais->Draw("triggeredPhiDeg-zoomPhiDeg>>"+name, "heading > -50 && l3TrigPatternH > 0", "goff");
    // h1->SetLineColor(dirInd+1);
    // h1->SetMaximum(600);
    // TString opt = dirInd == 0 ? "" : "sames";
    // h1->Draw(opt);

    TString name2 = TString::Format("h2_%d", dirInd);
    TH1D* h2 = new TH1D(name2,
    			varName + " " + titles[dirInd], 128, -5, 5);    
    // cWais->Draw("triggeredPhiDeg-phiExpected>>"+name2,
    cWais->Draw(varName + ">>"+name2,		
    		cuts, "goff");
    h2->SetLineColor(dirInd+1);
    TString opt = dirInd == 0 ? "" : "sames";
    h2->Draw(opt);

    h2->Scale(1./h2->Integral());

    hists.push_back(h2);
    
    // TString name3 = TString::Format("h3_%d", dirInd);
    // TH1D* h3 = new TH1D(name3,
    // 			varName + " " + titles[dirInd], 128, -5, 5);
    // cWais->Draw(varName + ">>"+name3, cuts, "goff");
    // if(dirInd==4){
    //   h3->SetLineColor(dirInd+2);
    // }
    // else{
    //   h3->SetLineColor(dirInd+1);
    // }
    // TString opt = dirInd == 0 ? "" : "sames";
    // // h3->Draw(opt);

    // hists.push_back(h3);

  }

  TCanvas* c1 = RootTools::drawHistsWithStatsBoxes(numDirs, &hists[0],
						   "", "mre");

  c1->BuildLegend();
  hists[0]->SetTitle("Debugging - run 352; #delta#phi(Degrees); Fraction of events/bin");
  hists[0]->SetMaximum(0.08);
}
