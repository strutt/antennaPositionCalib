void drawInputSummaryFile(){

  TString fName ="inputSummaryFile.root"; // single r, phi (grad, rms), t (mean), z(mean, rms, grad)

  TFile* f = TFile::Open(fName);

  TH1D* hPhi[NUM_SEAVEYS];
  // for(int ant=0; ant < NUM_SEAVEYS; ant++){  
  for(int phi=0; phi < NUM_PHI; phi++){
    new TCanvas();
    for(int ring=0; ring < 3; ring++){
      int ant = ring * NUM_PHI + phi;
      
      hPhi[ant] = (TH1D*) f->Get(TString::Format("hPhi%d", ant));

      hPhi[ant]->SetLineColor(ring + 1);
      TString opt = ring == 0 ? "" : "same";
      hPhi[ant]->Draw(opt);
    }
  }
}
