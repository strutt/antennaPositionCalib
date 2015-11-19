{

  const int numStages = 2;
  // TString stageNames[numStages] = {"before", "deltaT"};
  TString stageNames[numStages] = {"before", "deltaR"};  

  TString fileName = TString::Format("compareMeasuredWithExpectedDeltaTs_%sPlots.root", stageNames[1].Data());
  TFile* f = TFile::Open(fileName);

  TH1D h("h", "#deltat_{e} - #deltat_{m}; #phi_{expected} (Degrees); #deltat_{e} - #deltat_{m} (ns)", 360, 0, 360);

  const Int_t numColors = 4;

  auto cc = new CrossCorrelator();

  TCanvas* c1[numStages];
  TLegend* l[numStages];
  for(int stage=0; stage<numStages; stage++){
    c1[stage] = new TCanvas();
    l[stage] = new TLegend(0.87, 0, 1, 1);
    colCouter=0;
    for(int combo=0; combo<336; combo++){
      TString name = TString::Format("gr_%d_%s", combo, stageNames[stage].Data());
      TGraph* gr = (TGraph*) f->Get(name);

      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      if(ant1 < NUM_PHI && ant2 < NUM_PHI){
	if(gr){
	  switch (ant1%numColors){
	  case 0:
	    gr->SetMarkerColor(kRed);	    
	    break;
	  case 1:
	    gr->SetMarkerColor(kBlue);
	    break;
	  case 2:
	    gr->SetMarkerColor(kMagenta);
	    break;
	  case 3:
	    gr->SetMarkerColor(kCyan);
	    break;
	  }
	  if(ant1==0 && ant2==NUM_PHI-1){
	    gr->SetMarkerColor(kCyan);	    
	  }
	  gr->SetMarkerStyle(8);
	  l[stage]->AddEntry(gr, TString::Format("%d and %d", ant1, ant2), "p");
	  if(combo==0){
	    gr->Draw("ap");
	    gr->GetXaxis()->SetLimits(0, 360);
	  }
	  else{
	    gr->Draw("psame");
	  }
	}
      }
    }
    l[stage]->Draw();
  }
  

  

  
  return;

}
