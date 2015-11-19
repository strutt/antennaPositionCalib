{
  const Int_t run = 352;
  
  TString calEventFilePath = TString::Format("~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root",
					     run, run);


  // UInt_t eventNumber = 60832349; //Secondary peak event
  UInt_t eventNumber = 60832108; //Main peak event

  TFile* fE = TFile::Open(calEventFilePath);
  TTree* tE = (TTree*) fE->Get("eventTree");
  CalibratedAnitaEvent* calEvent = NULL;
  tE->SetBranchAddress("event", &calEvent);
  tE->BuildIndex("eventNumber"); 
  tE->GetEntryWithIndex(eventNumber);
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(calEvent);

  Bool_t savePlots = kFALSE;


  // TLine *l1 = new TLine(phiExpected, -75,
  // 			phiExpected,75);
  // TLine *l2 = new TLine(-11.25,thetaExpected,
  // 			360-11.25,thetaExpected);

  // l1->SetLineStyle(2);
  // l2->SetLineStyle(2);
  // TLegend* l0 = new TLegend(0.6, 0.91, 0.99, 0.99);
  // l0->AddEntry(l1, "WAIS divide location", "l");   
  
  
  // CrossCorrelator cc(32);
  CrossCorrelator cc;
  // cc.getNormalizedInterpolatedTGraphs(realEvent);
  cc.correlateEvent(realEvent);

  delete realEvent;


  TCanvas* c7 = new TCanvas();
  TH2D* hImageH = cc->makeImage(AnitaPol::kHorizontal);
  hImageH->Draw("colz");
  // l1->Draw();
  // l2->Draw();  
  // l0->Draw();
  
  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/imageH_event%u_WaisMarker", eventNumber);
    RootTools::saveCanvas(c7, saveName);
  }

  TCanvas* c8 = new TCanvas();
  TH2D* hImageV = cc->makeImage(AnitaPol::kVertical);  
  hImageV->Draw("colz");
  // l1->Draw();
  // l2->Draw();  
  // l0->Draw();

  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/imageV_event%u_WaisMarker", eventNumber);
    RootTools::saveCanvas(c8, saveName);
  }
}
