{


  TChain* t = new TChain("deltaTTree");
  const Int_t firstRun = 331;
  const Int_t lastRun = 354;
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("generateDeltaTTree_run%d-%dPlots.root", run, run);
    t->Add(fileName);
  }
  
  // TFile* f = TFile::Open("generateDeltaTTreePlots.root");
  // TTree* t = (TTree*) f->Get("deltaTTree");
  Bool_t savePlots = kFALSE;
  // Bool_t savePlots = kTRUE;  
  // return;

  // UInt_t eventNumber = 60832349; //Secondary peak event
  UInt_t eventNumber = 60832108; //Main peak event

  t->BuildIndex("eventNumber");
  Double_t thetaExpected;
  Double_t phiExpected;
  t->SetBranchAddress("thetaExpected", &thetaExpected);
  t->SetBranchAddress("phiExpected", &phiExpected);  
  t->GetEntryWithIndex(eventNumber);
  TGraph* grWaisExpected = new TGraph();
  phiExpected = phiExpected > 360 - 11.25 ? phiExpected - 360 : phiExpected;


  const Double_t deltaPhiMax = 45;
  const Double_t deltaPhiMin = -deltaPhiMax;

  TH2D* hDeltaTsVsPhi = new TH2D("hDeltaTsVsPhi",
  				 "#Deltat of maximum correlation between antennas 0 and 16 vs #delta#phi",
  				 // 100, deltaPhiMin, deltaPhiMax, 1024, -75, -75);
  				 1024, 0, 360, 100, -8, -3);
  hDeltaTsVsPhi->GetYaxis()->SetTitle("Correlation #Deltat (ns)");
  hDeltaTsVsPhi->GetXaxis()->SetTitle("#phi_{expected} - #phi_{phi-sector} (degrees)");

  // TCanvas* c3 = new TCanvas();
  // t->Draw("thetaExpected:phiExpected", //:correlationDeltaTs[0][6]",
  // 	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
  // 	  "colz");
  // t->Draw("thetaExpected:phiExpected", //:correlationDeltaTs[0][6]",
  // 	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
  // 	  "colz");
  t->Draw("thetaExpected:deltaPhiDeg[0]", //:correlationDeltaTs[0][6]",
	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
	  // "",
	  "colz");  

  // hDeltaTsVsPhi->Draw("colz");
  // hDeltaTsVsPhi->GetYaxis()->SetRangeUser(1, 9);
  // // c3->SetLogz(1);
  // if(savePlots){
  //   TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_vs_deltaPhi_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
  //   RootTools::saveCanvas(c3, fileName);
  // }


  return;
  
}
