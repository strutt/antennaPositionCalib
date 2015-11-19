{

  auto c = new TChain("angResTree");
  c->Add("hpolLindasHpolValues/generateAngularResolutionTreePlots_*.root");

  UInt_t l3TrigPatternH;
  Double_t heading;
  UInt_t eventNumber;  
  c->SetBranchAddress("heading", &heading);
  c->SetBranchAddress("l3TrigPatternH", &l3TrigPatternH);
  c->SetBranchAddress("eventNumber", &eventNumber);
  

  Long64_t numEntries = c->GetEntries();

  c->Show(numEntries-1);

  const UInt_t firstEvent = 55336937;
  const UInt_t lastEvent = 61535221;

  const Int_t numPhiDegBins = 64;
  const Int_t numEventNumberBins = 64;
  
  auto h = new TH2D("h", "Heading vs. L3 Triggered #phi-sectors for WAIS pulses; Heading (degrees); Triggered #phi-sector; Number of events", numPhiDegBins, 0, 360, NUM_PHI, 0, 16);
  auto h2 = new TH1D("h2", "Difference between heading and L3 Triggered #phi-sector position for WAIS pulses; Heading - Triggered #phi-sector position (Degrees); Number of events", numPhiDegBins, 0, 360);
  auto h3 = new TH2D("h3", "Difference between heading and L3 Triggered #phi-sector position for WAIS pulses; Event Number; Heading - Triggered #phi-sector position (Degrees); Number of events", numEventNumberBins, firstEvent, lastEvent+1, numPhiDegBins, 0, 360);

  auto geom = AnitaGeomTool::Instance();
  
  for(Long64_t entry=0; entry < numEntries; entry++){
    c->GetEntry(entry);

    if(heading > -20){
      int numTrig=0;
      for(int phi=0; phi<NUM_PHI; phi++){
	if(RootTools::getBit(phi, l3TrigPatternH)){

	  Double_t antPhiDeg = geom->getAntPhiPositionRelToAftFore(phi)*TMath::RadToDeg();

	  h->Fill(heading, phi);
	  
	  Double_t dPhiDeg = RootTools::getDeltaAngleDeg(heading, antPhiDeg);

	  dPhiDeg = dPhiDeg < 0 ? dPhiDeg + 360 : dPhiDeg;

	  h2->Fill(dPhiDeg);
	  h3->Fill(eventNumber, dPhiDeg);
	  numTrig++;
	}
      }
      // if(numTrig==0){
      // 	std::cout << entry << "\t" << heading << "\t" << l3TrigPatternH << std::endl;
      // }
    }
  }
  
  auto c1 = new TCanvas();
  h->Draw("colz");
  h->GetXaxis()->SetNoExponent(1);
  h->GetYaxis()->SetNoExponent(1);      
  h->GetZaxis()->SetNoExponent(1);
  auto c2 = new TCanvas();
  h2->Draw();
  h2->GetXaxis()->SetNoExponent(1);
  h2->GetYaxis()->SetNoExponent(1);
  h2->GetZaxis()->SetNoExponent(1);      
  auto c3 = new TCanvas();
  h3->Draw("colz");
  h3->GetXaxis()->SetNoExponent(1);
  h3->GetYaxis()->SetNoExponent(1);
  h3->GetZaxis()->SetNoExponent(1);    









}
