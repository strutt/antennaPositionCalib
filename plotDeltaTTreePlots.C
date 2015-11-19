plotDeltaTTreePlots(){

  gSystem->Load("libBensAnitaTools.so");
  
  TFile* f = TFile::Open("generateDeltaTTreePlots.root");
  TTree* t = (TTree*) f->Get("deltaTTree");

  CrossCorrelator* cc = new CrossCorrelator();


  Int_t combo=6;
  for(combo=combo; combo < 366; combo++){
    Int_t ant1 = cc->comboToAnt1s.at(combo);
    Int_t ant2 = cc->comboToAnt2s.at(combo);
    if(TMath::Abs(RootTools::getDeltaAngleDeg(cc->phiArrayDeg[ant1], cc->phiArrayDeg[ant2])) < 5){
      break;
    }
  }
  cout << combo << "\t" << ant1 << "\t" << ant2 << std::endl;
  cout << combo << "\t" << cc->phiArrayDeg[ant1] << "\t" << cc->phiArrayDeg[ant2] << endl;
  

  // TString command = TString::Format("correlationDeltaTs[0][%d]:phiExpected", combo);
  TString command = TString::Format("phiExpected", combo);  
  // TString command = TString::Format("thetaExpected:phiExpected", combo);  
  TString opt = "";
  // opt += TString::Format("phiExpected > 315 - 45 && phiExpected < 315 + 45 "); //, ant1);
  // opt += TString::Format("&& TMath::Abs(deltaPhiDeg[%d]) < 33.375", ant2);
  // t->Draw(command, "", "colz"); //, opt);
  t->Draw(command, opt, "colz");
  
  delete cc;
  
}
