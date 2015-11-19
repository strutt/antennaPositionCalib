void binSoftwareTriggersByL3TriggersH(){

  auto c = RootTools::getHeadChain(127, 439);

  RawAnitaHeader* header = nullptr;
  c->SetBranchAddress("header", &header);

  Long64_t nEntries = c->GetEntries();
  std::cout << nEntries << std::endl;
  const Int_t numPhi = 16;
  ProgressBar p(nEntries);

  const UInt_t firstEvent = 9242201;
  const UInt_t lastEvent = 10809498;
  const Int_t numEventBins = 512;
  
  auto hL3Triggers = new TH1D("hL3Triggers", "Distribution of L3 VPOL triggered #Phi-sectors (trigType & 1)==0 (Non-RF triggers); #Phi-sector; Number of events", numPhi, 0, numPhi);
  hL3Triggers->GetYaxis()->SetNoExponent(1);
  auto hL3Triggers2 = new TH2D("hL3TriggersH", "Distribution of L3 VPOL triggered #Phi-sectors (trigType & 1)==0 (Non-RF triggers); eventNumber; #Phi-sector; Number of events", numEventBins, firstEvent, lastEvent+1, numPhi, 0, numPhi);
  hL3Triggers2->GetXaxis()->SetNoExponent(1);

  for(Long64_t entry=0; entry<nEntries; entry++){
    
    c->GetEntry(entry);

    if((header->trigType & 1)==0){ // Software trigger

      if(entry==0 || entry == nEntries - 1){
	c->Show(entry);
      }

      for(UInt_t phiSect=0; phiSect<numPhi; phiSect++){
	Int_t triggered = RootTools::getBit(phiSect, header->l3TrigPatternH);
	if(triggered > 0){
	  hL3Triggers->Fill(phiSect);
	  hL3Triggers2->Fill(header->eventNumber, phiSect);
	}
      }
    }
    p++;
  }

  auto c1 = new TCanvas();
 
  hL3Triggers->Draw();
  hL3Triggers->SetMinimum(0);  
  
  auto c2 = new TCanvas();
  hL3Triggers2->Draw("colz");  
}
