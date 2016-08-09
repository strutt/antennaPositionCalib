void draw5(){

  gStyle->SetOptStat("mre");


  // TFile* f = TFile::Open("testing/corInterp/latestNums/generateAngularResolutionTreePlots_352_2016-02-21_15-38-05.root"); //generateAngularResolutionTreePlots_352_2016-02-19_18-01-52.root");
  TChain* c = new TChain("angResTree");

  UInt_t l3TrigPatternH;
  c->SetBranchAddress("l3TrigPatternH", &l3TrigPatternH);

  Double_t phiExpected;
  c->SetBranchAddress("phiExpected", &phiExpected);  
  // c->Add("newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15/*.root");
  c->Add("newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15/generateAngularResolutionTreePlots_352*.root");  

  const Long64_t numEntries = c->GetEntries();
  ProgressBar p(numEntries);

  TH2D* h = new TH2D("h", "#phi_{expected} vs. L3 triggered #phi-sectors; #phi_{expected} (Degrees); #phi-sector", 360, 0, 360, NUM_PHI, 0, NUM_PHI);
  TH1D* h2 = new TH1D("h2", "Angle between triggered #phi-sector center and #phi_{expected}; #delta#phi (Degrees); Events / bin", 360, -180, 180);
  TH1D* h3 = new TH1D("h3", "#phi-sectors triggered; #phi-sector; Number of triggers", NUM_PHI, 0, NUM_PHI);
  TH1D* h4 = new TH1D("h4", "#delta#phi-sectors; #phi-sector; Number of triggers", 5, -2, 3);

  for(Long64_t entry=0; entry < numEntries; entry++){

    c->GetEntry(entry);

    for(int phi=0; phi<NUM_PHI; phi++){
      Int_t l3 = RootTools::getBit(phi, l3TrigPatternH);
      if(l3){
	h->Fill(phiExpected, phi);

	Double_t deltaPhiThingy = RootTools::getDeltaAngleDeg(phiExpected + 45, 22.5*phi);
	h2->Fill(deltaPhiThingy);

	h3->Fill(phi);

	if(deltaPhiThingy <= -11.25){
	  h4->Fill(-1);
	}
	else if(deltaPhiThingy > 11.25){
	  h4->Fill(1);
	}
	else{
	  h4->Fill(0);
	}
      }
    }
    p++;
  }
  new TCanvas();
  h->Draw("colz");
  new TCanvas();
  h2->Draw();
  new TCanvas();
  h4->Draw();
  for(int binx=1; binx<=h4->GetNbinsX(); binx++){
    h4->GetXaxis()->SetBinLabel(binx, TString::Format("%d", -3 + binx));
  }
  h4->GetXaxis()->SetLabelOffset(.01);
    
}
