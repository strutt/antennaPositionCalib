void draw3(){

  gStyle->SetOptStat("mre");
  
  const int numDirs = 2;
  // TString dirNames[numDirs]={"testing/corInterp",
  TString dirNames[numDirs]={"newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15",			     
			     "../cosminsTrees"};

  // TString titles[numDirs] = {"bf1", "bins0.5", "bigBins", "corInterp", "corInterp_photo"};
  
  TString titles[numDirs] = {"Ben's reco", "Cosmin's reco"};


  std::vector<TH1D*> hists;


  TChain* cs[numDirs];

   // auto c1 = new TCanvas();
  for(int dirInd=0; dirInd < numDirs; dirInd++){
    auto cWais = new TChain("angResTree");
    cWais->Add(dirNames[dirInd] + "/generateAngularResolutionTreePlots_352*root");

    
    TString varName = "deltaPhiDeg";
    TString cuts = "heading > -50 && l3TrigPatternH > 0";
    
    if(cWais->GetEntries() == 0){
      cWais = new TChain("pointed");
      cWais->Add(dirNames[dirInd] + "/352.root");      
      varName = "-1*WAIS.dphi";
      cuts = "heading > -50";

    }

    cWais->BuildIndex("eventNumber");
    cs[dirInd] = cWais;
    
    const double angleRange = 5; //3; //180;
    const int numBins = 64;

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

    // TString name2 = TString::Format("h2_%d", dirInd);
    // TH1D* h2 = new TH1D(name2,
    // 			"triggeredPhiDeg-phiExpected: " + titles[dirInd], 128, -5, 5);    
    // cWais->Draw("triggeredPhiDeg-phiExpected>>"+name2,
    // 		"heading > -50 && l3TrigPatternH > 0", "goff");
    // h2->SetLineColor(dirInd+1);
    // TString opt = dirInd == 0 ? "" : "sames";
    // h2->Draw(opt);
    
    TString name3 = TString::Format("h3_%d", dirInd);
    TH1D* h3 = new TH1D(name3,
    			// varName + " " + titles[dirInd], numBins, -angleRange, angleRange);
    			titles[dirInd], numBins, -angleRange, angleRange);    
    cWais->Draw(varName + ">>"+name3, cuts, "goff");
    if(dirInd==4){
      h3->SetLineColor(dirInd+2);
    }
    else{
      h3->SetLineColor(dirInd+1);
    }
    TString opt = dirInd == 0 ? "" : "sames";
    // h3->Draw(opt);

    hists.push_back(h3);

  }

  TCanvas* c1 = RootTools::drawHistsWithStatsBoxes(numDirs, &hists[0],
						   "", "mre");

  new TCanvas();
  c1->BuildLegend();
  hists[0]->SetTitle("Debugging - run 352; #delta#phi(Degrees); Number of events");

  cs[0]->AddFriend(cs[1]);
  
  const int numEntries = cs[0]->GetEntries();

  TH2D* hTest = new TH2D("h2", "h2", 64, -3, 3, 64, -3, 3);
  // TH2D* hTest = new TH2D("h2", "h2", 256, 0, 360, 256, 0, 360);  

  // cs[0]->Draw("deltaPhiDeg:WAIS.dphi>>h2", "", "colz");
  // cs[0]->Draw("deltaPhiDeg:-WAIS.dphi>>h2", "", "colz");
  cs[0]->Draw("deltaPhiDeg:-1*WAIS.dphi>>h2", "TMath::Abs(deltaThetaDeg) < 5", "colz");  
  hTest->SetTitle("Ben New vs. Cosmin (run 352); -1#timesWAIS.dphi (Degrees); deltaPhiDeg (Degrees)");
  
  // cs[0]->Draw("deltaPhiDeg + WAIS.dphi", "TMath::Abs(WAIS.dphi) < 5 && TMath::Abs(deltaPhiDeg) < 5 ", "colz");
  // cs[0]->Draw("deltaPhiDeg - WAIS.dphi", "TMath::Abs(WAIS.dphi) < 5 && TMath::Abs(deltaPhiDeg) < 5 ", "colz");    

  // cs[0]->Draw("phiExpected:WAIS.phi>>h2", "", "colz");

  // cs[0]->Draw("phiExpected:WAIS.phi>>h2", "", "colz");
  new TCanvas();
  // cs[0]->Draw("phiExpected-WAIS.phi", "", "");

  auto h333 = new TH1D("h333", "Cosmin's vs. Ben's #phi_{expected}; #delta#phi_{expected} (Degrees); Number of events / bin", 128, -5, 5);
  cs[0]->Draw("phiExpected-WAIS.phi>>h333", "", "");
  
  
  // for(int entry=0; entry < numEntries; entry++){

  //   cs[0]->GetEntry(entry);
  //   // std::cout << cs[1]->GetEntryWithIndex(eventNumber) << std::endl;
  //   // cs[1]->GetEntryWithIndex(eventNumber);
  //   cs[1]->GetEntry(entry);

  //   // std::cout << deltaPhiDeg << "\t" << dphi << std::endl;
    
  //   hTest->Fill(deltaPhiDeg, dphi);
    
  // }
  // new TCanvas();
  // hTest->Draw("colz");

}
