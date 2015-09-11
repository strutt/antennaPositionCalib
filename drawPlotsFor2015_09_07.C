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
  
  Double_t length = 5;
  TLine *l1 = new TLine(phiExpected, -75,
  			phiExpected,75);
  TLine *l2 = new TLine(-11.25,thetaExpected,
  			360-11.25,thetaExpected);

  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  TLegend* l0 = new TLegend(0.6, 0.91, 0.99, 0.99);
  l0->AddEntry(l1, "WAIS divide location", "l");   
  
  TH1D* hDeltaTs = new TH1D("hDeltaTs",
			    "#Deltat of maximum correlation between antennas 0 and 16",
  			    1024, -75, -75);
  hDeltaTs->GetXaxis()->SetTitle("Correlation #Deltat (ns)");
  hDeltaTs->GetYaxis()->SetTitle("Number of events");  
  // t->Draw("correlationDeltaTs[0][6]>>hDeltaTs", "TMath::Abs(deltaPhiDeg[0])<22.5", "goff");
  t->Draw("correlationDeltaTs[0][6]>>hDeltaTs",
	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
	  "goff");  
  TCanvas* c1 = new TCanvas();
  hDeltaTs->Draw();
  c1->SetLogy(1);
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_waisSelection_runs%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c1, fileName);
  }

  
  TH1D* hDeltaTsZoom = hDeltaTs->Clone("hDeltaTsZoom");
  TCanvas* c2 = new TCanvas();
  hDeltaTsZoom->Draw();
  hDeltaTsZoom->GetXaxis()->SetRangeUser(-3, 12);
  c2->SetLogy(1);
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_zoom_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c2, fileName);
  }
  
  return;

  TH2D* hDeltaTsVsCorr = new TH2D("hDeltaTsVsCorr",
				  "#Deltat of maximum correlation between antennas 0 and 16 vs correlation",
				  128, 0, 1, 1024, -75, -75);
  hDeltaTsVsCorr->GetYaxis()->SetTitle("Correlation #Deltat (ns)");
  hDeltaTsVsCorr->GetXaxis()->SetTitle("Correlation at peak #deltat (no units)");

  TCanvas* c2a = new TCanvas();
  t->Draw("correlationDeltaTs[0][6]:correlationValues[0][6]>>hDeltaTsVsCorr",
	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
	  "goff");
  hDeltaTsVsCorr->Draw("colz");

  TH2D* hDeltaTsVsCorrZoom = (TH2D*) hDeltaTsVsCorr->Clone();
  hDeltaTsVsCorr->RebinY(8);
  
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_vs_corrVal_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c2a, fileName);
  }
  TCanvas* c2b = new TCanvas();
  hDeltaTsVsCorrZoom->GetYaxis()->SetRangeUser(1, 9);
  hDeltaTsVsCorrZoom->Draw("colz");
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_vs_corrVal_zoom_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c2b, fileName);
  }




  
  

  TH2D* hDeltaTsVsPhi = new TH2D("hDeltaTsVsPhi",
			    "#Deltat of maximum correlation between antennas 0 and 16 vs #delta#phi",
  			    100, deltaPhiMin, deltaPhiMax, 1024, -75, -75);
  hDeltaTsVsPhi->GetYaxis()->SetTitle("Correlation #Deltat (ns)");
  hDeltaTsVsPhi->GetXaxis()->SetTitle("#phi_{expected} - #phi_{phi-sector} (degrees)");

  TCanvas* c3 = new TCanvas();
  t->Draw("correlationDeltaTs[0][6]:deltaPhiDeg[0]>>hDeltaTsVsPhi",
	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
	  "goff");
  hDeltaTsVsPhi->Draw("colz");
  hDeltaTsVsPhi->GetYaxis()->SetRangeUser(1, 9);
  // c3->SetLogz(1);
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_vs_deltaPhi_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c3, fileName);
  }


  TH2D* hDeltaTsVsTheta = new TH2D("hDeltaTsVsTheta",
			    "#Deltat of maximum correlation between antennas 0 and 16 vs #theta_{expected}",
  			    100, deltaPhiMin, deltaPhiMax, 1024, -75, -75);
  hDeltaTsVsTheta->GetYaxis()->SetTitle("Correlation #Deltat (ns)");
  hDeltaTsVsTheta->GetXaxis()->SetTitle("#theta_{expected} (degrees)");

  TCanvas* c4 = new TCanvas();
  t->Draw("correlationDeltaTs[0][6]:thetaExpected>>hDeltaTsVsTheta",
	  TString::Format("TMath::Abs(deltaPhiDeg[0])<%lf", deltaPhiMax),
	  "goff");
  hDeltaTsVsTheta->Draw("colz");
  hDeltaTsVsTheta->GetYaxis()->SetRangeUser(1, 9);
  // c3->SetLogz(1);
  if(savePlots){
    TString fileName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/maxCorrDt_vs_thetaExpected_waisSelection_run%d-%d_ants_0_16", firstRun, lastRun);
    RootTools::saveCanvas(c4, fileName);
  }


  return;
  

  TString calEventFilePath = "~/UCL/ANITA/calibratedFlight1415/run352/calEventFile352.root";
  TFile* fE = TFile::Open(calEventFilePath);
  TTree* tE = (TTree*) fE->Get("eventTree");
  CalibratedAnitaEvent* calEvent = NULL;
  tE->SetBranchAddress("event", &calEvent);
  tE->BuildIndex("eventNumber"); 
  tE->GetEntryWithIndex(eventNumber);
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(calEvent);

  CrossCorrelator cc(32);
  // CrossCorrelator cc(1);  
  // cc.getNormalizedInterpolatedTGraphs(realEvent);
  cc.correlateEvent(realEvent);

  delete realEvent;

  
  
  TCanvas* c5 = new TCanvas();
  const Int_t numGrs = 2;

  TGraph* grs[numGrs];
  TString legTitles[numGrs];
  Int_t ants[numGrs] = {0, 16};
  Double_t rmsVals[numGrs];
  TString grTitle = TString::Format("Interpolated TGraphs for event %u for antennas",eventNumber);
  for(Int_t antInd=0; antInd < numGrs; antInd++){
    // grs[antInd] = realEvent->getGraph(ants[antInd], AnitaPol::kHorizontal);
    grs[antInd] = (TGraph*) cc->grsInterp[AnitaPol::kHorizontal][ants[antInd]]->Clone();

    rmsVals[antInd] = cc->interpRMS[AnitaPol::kHorizontal][ants[antInd]];
    legTitles[antInd] = TString::Format("Antenna %d", ants[antInd]);
    grTitle += TString::Format(" %d", ants[antInd]);
    if(antInd != numGrs-1){
      grTitle += TString::Format(",");
    }
  }
  RootTools::multiplyTGraphYAxes(numGrs, grs, rmsVals);
  grs[0]->SetTitle(grTitle);
  grs[0]->GetXaxis()->SetTitle("Time (ns)");
  grs[0]->GetXaxis()->SetRangeUser(-4, 100);
  grs[0]->GetYaxis()->SetTitle("Voltage (mV)");
  
  
  RootTools::drawArrayOfTGraphsPrettily(grs, numGrs, "l", c5, NULL);
  TLegend* l = RootTools::makeLegend(grs, numGrs, legTitles, "l", 0.8, 0.8, 1, 1);
  l->Draw();
  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/interp_TGraphs_waisSelection_secondary_peak_run352_event%u_ants_0_16", eventNumber);
    RootTools::saveCanvas(c5, saveName);
  }

  TCanvas* c6 = new TCanvas();
  TGraph* grsOffset[numGrs];
  Double_t offsets[numGrs] = {0, 7.4399038};
  grTitle = TString::Format("Correlation offset interpolated TGraphs for event %u for antennas",eventNumber);
  for(Int_t antInd=0; antInd < numGrs; antInd++){
    grsOffset[antInd] = (TGraph*) grs[antInd]->Clone();
    grTitle += TString::Format(" %d", ants[antInd]);
    if(antInd != numGrs-1){
      grTitle += TString::Format(",");
    }

  }
  RootTools::offsetTGraphXAxes(numGrs, grsOffset, offsets);
  RootTools::drawArrayOfTGraphsPrettily(grsOffset, numGrs, "l", c6, NULL);
  l->Draw();
  grsOffset[0]->SetTitle(grTitle);

  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/interp_TGraphs_offset_waisSelection_secondary_peak_run352_event%u_ants_0_16", eventNumber);
    RootTools::saveCanvas(c6, saveName);
  }


  TCanvas* c7 = new TCanvas();
  TH2D* hImageH = cc->makeImage(AnitaPol::kHorizontal);
  hImageH->Draw("colz");
  l1->Draw();
  l2->Draw();  
  l0->Draw();
  

  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/imageH_event%u_WaisMarker", eventNumber);
    RootTools::saveCanvas(c7, saveName);
  }

  TCanvas* c8 = new TCanvas();
  TH2D* hImageV = cc->makeImage(AnitaPol::kVertical);  
  hImageV->Draw("colz");
  l1->Draw();
  l2->Draw();  
  l0->Draw();

  if(savePlots){
    TString saveName = TString::Format("~/UCL/ANITA/weeklyMeeting/2015_09_07/imageV_event%u_WaisMarker", eventNumber);
    RootTools::saveCanvas(c8, saveName);
  }

}
