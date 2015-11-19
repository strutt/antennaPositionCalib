void binWaisPulsesByL3TriggersH(){

  auto c = new TChain("angResTree");
  c->Add("hpolLindasHpolValues/*.root");

  UInt_t l3TrigPatternH;
  c->SetBranchAddress("l3TrigPatternH", &l3TrigPatternH);
  UInt_t eventNumber;
  c->SetBranchAddress("eventNumber", &eventNumber);
  Double_t deltaPhiDeg;
  c->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);

  Long64_t nEntries = c->GetEntries();
  std::cout << nEntries << std::endl;
  const Int_t numPhi = 16;
  ProgressBar p(nEntries);

  const UInt_t firstEvent = 55336937;
  const UInt_t lastEvent = 61535221;
  const Int_t numEventBins = 512;
  
  auto hL3TriggersH = new TH1D("hL3TriggersH", "Distribution of L3 triggered #Phi-sectors (WAIS pulses); #Phi-sector; Number of events", numPhi, 0, numPhi);
  hL3TriggersH->GetYaxis()->SetNoExponent(1);
  auto hL3TriggersH2 = new TH2D("hL3TriggersH", "Distribution of L3 triggered #Phi-sectors (WAIS pulses); eventNumber; #Phi-sector; Number of events", numEventBins, firstEvent, lastEvent+1, numPhi, 0, numPhi);
  hL3TriggersH2->GetXaxis()->SetNoExponent(1);
  
  for(Long64_t entry=0; entry<nEntries; entry++){
    
    c->GetEntry(entry);

    if(entry==0 || entry == nEntries - 1){
      c->Show(entry);
    }

    for(Int_t phiSect=0; phiSect<numPhi; phiSect++){
      Int_t triggered = RootTools::getBit(phiSect, l3TrigPatternH);

      if(triggered > 0 && TMath::Abs(deltaPhiDeg) < 5){
	hL3TriggersH->Fill(phiSect);
	hL3TriggersH2->Fill(eventNumber, phiSect);
      }
    }
    p++;
  }

  auto c1 = new TCanvas();
 
  hL3TriggersH->Draw();
  hL3TriggersH->SetMinimum(0);  
  
  auto c2 = new TCanvas();
  hL3TriggersH2->Draw("colz");  
}
