void drawMakeCorrelationSummariesTimePlots(){

  gStyle->SetOptStat("mre");

  // TFile* f = TFile::Open("makeCorrelationSummariesTimePlots_352_2016-02-26_18-02-20.root");
  auto corTree = new TChain("corTree");
  corTree->Add("makeCorrelationSummariesTimePlots/*.root");
  // TTree* corTree = (TTree*) f->Get("corTree");
  cout << corTree << endl;

  auto cc = new CrossCorrelator();
  auto h = new TH1D("h", "", 128, -2, 2);

  // Book histogram and set x, y axis labels to antenna names
  TString name = "h2";
  TString polLetter = "H";
  // TString polLetter = pol == AnitaPol::kHorizontal ? "H" : "V";  
  // name += polLetter;
  // name += TString::Format("%u", eventNumber[pol]);

  TString title = TString::Format("Mean #deltat_{measured} - #deltat_{expected} for all used antenna pairs; Antenna 1; Antenna 2; Mean #Delta#deltat (ns)");
  
  int binShift = 0;
  // title += pol==AnitaPol::kHorizontal ? "HPOL" : "VPOL";
  // title += "; ; ; Correlation coefficient";
  TH2D* h2 = new TH2D(name, title, NUM_SEAVEYS, 0, NUM_SEAVEYS, NUM_SEAVEYS, 0, NUM_SEAVEYS);
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    TString lab = TString::Format("%d", 1 + (ant/NUM_RING));
    if(ant%NUM_RING==0){
      lab += "T"+polLetter;
    }
    else if((ant%NUM_RING)==1){
      lab += "M"+polLetter;
    }
    else{
      lab += "B"+polLetter;
    }
    int binInd = ant + binShift;
    // binInd = binInd < 0 ? binInd + NUM_SEAVEYS : binInd;
    binInd = binInd >= NUM_SEAVEYS ? binInd - NUM_SEAVEYS : binInd;    
    h2->GetXaxis()->SetBinLabel(binInd+1, lab.Data());
    h2->GetYaxis()->SetBinLabel(binInd+1, lab.Data());
  }

  std::vector<TH1D*> hs;
  for(int combo=0; combo < NUM_COMBOS; combo++){

    int ant1 = cc->comboToAnt1s.at(combo);
    int ant2 = cc->comboToAnt2s.at(combo);

    auto command = TString::Format("deltaTMeasured[%d]-deltaTExpected[%d]>>h", combo, combo);
    auto cuts = TString::Format("deltaTMeasured[%d]>-999 && TMath::Abs(deltaTMeasured[%d]-deltaTExpected[%d]) < 0.7", combo, combo, combo);
    // cout << command << "\t" << cuts << endl;

    corTree->Draw(command, cuts, "goff");
    
    std::cout << ant1 << "\t" << ant2<< "\t" << h->GetMean() << std::endl;

    // auto h3 = ((TH1D*)h->Clone());
    // h3->SetTitle(TString::Format("antennas %d and %d; #Delta#deltat (ns); Events / bin", ant1, ant2));

    // hs.push_back(h3);

    Int_t phi1 = ant1%NUM_PHI;
    Int_t phi2 = ant2%NUM_PHI;

    Int_t ring1 = ant1/NUM_PHI;
    Int_t ring2 = ant2/NUM_PHI;

    Int_t theBinX = phi1*NUM_RING + ring1 + 1;
    Int_t theBinY = phi2*NUM_RING + ring2 + 1;

    theBinX += binShift;
    theBinY += binShift;

    theBinX = theBinX >= NUM_SEAVEYS ? theBinX - NUM_SEAVEYS : theBinX;
    theBinY = theBinY >= NUM_SEAVEYS ? theBinY - NUM_SEAVEYS : theBinY;
      
    h2->SetBinContent(theBinX, theBinY, h->GetMean());
    h2->SetBinContent(theBinY, theBinX, h->GetMean());          
    
    // if(combo >= 5){
    //   break;
    // }
  }
  
  // TCanvas* c1 = RootTools::drawArrayOfHistosPrettily(&hs[0], (Int_t) hs.size());
  // TCanvas* c2 = RootTools::drawHistsWithStatsBoxes((Int_t) hs.size(), &hs[0],
  // 						   "", "mre");
  // TLegend* l = c2->BuildLegend();
  // hs[0]->SetTitle("#deltat_{measured} - #deltat_{expected} for a selection of antenna pairs; #Delta#deltat (ns); Events / bin");
  
  auto c3 = new TCanvas();
  h2->Draw("colz");
}
