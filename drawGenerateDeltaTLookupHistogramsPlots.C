void drawGenerateDeltaTLookupHistogramsPlots(){
 
  
  TFile* f = TFile::Open("generateDeltaTLookupHistogramsPlots.root");

  RootTools::setWhiteZeroColorScale();
  
  std::vector<Int_t> combos;
  std::vector<Int_t> ant1s;
  std::vector<Int_t> ant2s;  
  CrossCorrelator* cc = new CrossCorrelator();

  const Int_t numCombos = 336; //NUM_PHI;
  const Int_t numPhi = 16;

  // // Draw feed and photogrammetry examples for 0, 16.
  // TH2D* hPureDiff_2D_feed_0_16 = (TH2D*) f->Get("hPureDiff_2D_feed_0_16");
  // TH2D* hPureDiff_2D_photo_0_16 = (TH2D*) f->Get("hPureDiff_2D_photo_0_16");
  // RootTools::draw2D(hPureDiff_2D_feed_0_16, "colz");
  // return;

  // for(Int_t comboInd=0; comboInd < numCombos; comboInd++){
  // Int_t startCombo = 0;
  Int_t startCombo = 0; //NUM_COMBOS - 50;
  Int_t endCombo = 1; //NUM_COMBOS;

  startCombo = 0;
  endCombo = 50;
  // Int_t endCombo = 50;
  for(Int_t comboInd=startCombo; comboInd < endCombo; comboInd++){

    Int_t combo = comboInd;
    Int_t ant1 = cc->comboToAnt1s.at(comboInd);
    Int_t ant2 = cc->comboToAnt2s.at(comboInd);


    if(!(TMath::Abs(ant1 - ant2) == NUM_PHI || TMath::Abs(ant1 - ant2) == (2*NUM_PHI) || TMath::Abs(ant1 - ant2) == 1 || TMath::Abs(ant1 - ant2) == (NUM_PHI-1))){
      continue;
    }

    // int selPhi = 0;
    // if(!(ant1==selPhi || ant1==selPhi+NUM_PHI)){
    //   continue;
    // }

    // Int_t phiSect1 = ant1%numPhi;
    // Int_t phiSect2 = ant2%numPhi;
    // Int_t deltaPhiSect = phiSect1 - phiSect2;
    // cout << deltaPhiSect << "\t" << phiSect1 << "\t "<< phiSect2 << endl;
    // if(TMath::Abs(deltaPhiSect) > 1 && TMath::Abs(deltaPhiSect) < 15){
    //   continue;
    // }
    
    
    // Int_t ant1 = comboInd; // i.e. phi
    // // Int_t ant2 = (comboInd + 1)%16; //NUM_PHI;    
    // // Int_t ant2 = comboInd + 2*16; //NUM_PHI;
    // Int_t ant2 = comboInd + 16; //NUM_PHI;
    // Int_t combo = cc->comboIndices[ant1][ant2];
    // // std:: cout << ant1 << "\t" << ant2 << "\t" << combo << std::endl;
    combos.push_back(combo);
    ant1s.push_back(ant1);
    ant2s.push_back(ant2);

    TString name = TString::Format("hDtSparse_%d_%d", ant1, ant2);
    THnSparseF* hDtSparse = (THnSparseF*) f->Get(name);
    TCanvas* c1 = new TCanvas();
    c1->Divide(2);
    c1->cd(1);
    TH2D* hDtProf = hDtSparse->Projection(1, 0);
    RootTools::draw2D(hDtProf, "colz");
    hDtProf->GetYaxis()->SetRangeUser(-9, -5);
        

    c1->cd(2);
    name = TString::Format("hCorrDts_%d_%d", ant1, ant2);
    THnSparseF* hCorrDts = (THnSparseF*) f->Get(name);
    // TCanvas* c1 = new TCanvas();
    TH1D* hCorrDts_px = hCorrDts->Projection(0);
    hCorrDts_px->GetYaxis()->SetTitle("Number of events");
    hCorrDts_px->GetYaxis()->SetNoExponent(1);    
    hCorrDts_px->Draw();
    hCorrDts_px->GetXaxis()->SetRangeUser(-10, 10);
    gPad->SetLogy(1);

    cout << combo << "\t" << ant1 << "\t" << ant2 << "\t" << hCorrDts_px->GetMean() << endl;


    
  }

  delete cc;
  
}


