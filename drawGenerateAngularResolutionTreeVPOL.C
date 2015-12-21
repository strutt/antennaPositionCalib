#include "FFTtools.h"

void drawGenerateAngularResolutionTreeVPOL(){
  // gStyle->SetOptFit(0000);
  gStyle->SetOptStat(1111);  

  auto f = TFile::Open("newLindaNumbers_4steps_VPOL_10kVSeavey_2015_12_07_time_18_18_53/generateAngularResolutionTreeVPOLPlots_152_2015-12-19_00-08-43.root");
  auto c1 = new TCanvas();
  auto t = (TTree*) f->Get("angResTree");
  t->Draw("deltaPhiDeg");
  Double_t deltaPhiDeg;
  UInt_t eventNumber;  
  t->SetBranchAddress("deltaPhiDeg", &deltaPhiDeg);
  t->SetBranchAddress("eventNumber", &eventNumber);  
  

  auto fHeader = TFile::Open("~/UCL/ANITA/flight1415/root/run152/headFile152.root");
  auto tHeader = (TTree*) fHeader->Get("headTree");
  RawAnitaHeader* header = NULL;
  tHeader->BuildIndex("eventNumber");
  tHeader->SetBranchAddress("header", &header);
  
  auto fPat = TFile::Open("~/UCL/ANITA/flight1415/root/run152/gpsFile152.root");
  auto tPat = (TTree*) fPat->Get("adu5PatTree");
  tPat->BuildIndex("realTime");
  Adu5Pat* pat = NULL;
  tPat->SetBranchAddress("pat", &pat);

  // tPat->Show(0);
  // tPat->Show(tPat->GetEntries()-1);

  Double_t minTime = 1418940809;
  Double_t maxTime = 1418945220;

  const Int_t numTimeBins = 2048;

  // TH2D* h = new TH2D("h", "h", numTimeBins, -10, 10, 2, 0, 2);
  TH2D* h = new TH2D("h", "h", numTimeBins, minTime, maxTime+1, numTimeBins, -10, 10); 
  TH2D* hRawHeading = new TH2D("hRawHeading", "hRawHeading", numTimeBins, minTime, maxTime+1, numTimeBins, 0, numTimeBins);      

  Double_t lastHeading = 0;
  

  tPat->GetEntry(0);
  Double_t firstHeading = pat->heading;
  
  for(Long64_t entry=0; entry<tPat->GetEntries(); entry++){
    tPat->GetEntry(entry);

    

    Double_t deltaHeading = pat->heading - lastHeading;    

    if(deltaHeading > 360){
      deltaHeading -= 360;
    }
    else if(deltaHeading < -360){
      deltaHeading += 360;
    }

    h->Fill(pat->realTime, deltaHeading);
       
    
    lastHeading = pat->heading;
  }
  
  // tPat->Draw("pat->attFlag");


  h->Draw("colz");
  auto h2 = h->ProfileX();
  h2->Draw("same");


  
  // TF1* f1 = new TF1("fit", TString::Format("[0]+[1]*sin([2]*(x-%lf) + [3])", h->GetXaxis()->GetBinLowEdge(1));)
  TF1* f1 = new TF1("fit", TString::Format("[0]+[1]*sin([2]*(x-%lf) + [3]) + [4]*sin([5]*(x-%lf) + [6])",
					   minTime, minTime));

  // integral of a*sin(b*x + c) = -(a*sin(b*x+c))/b

  TF1* f1Int = new TF1("fitInt", TString::Format("[7] + [0]*(x-%lf) -[1]*cos([2]*(x-%lf)+[3])/[2] - [4]*cos([5]*(x-%lf)+[6])/[5]", minTime, minTime, minTime));

  
  f1->SetParameter(0, h->ProjectionY()->GetMean());
  f1->SetParameter(1, 0.2);
  // f1->SetParameter(2, 1./120);
  f1->SetParameter(2, TMath::TwoPi()/120);
  f1->SetParameter(3, TMath::Pi());
  f1->SetParameter(4, 0.1);
  f1->SetParameter(5, TMath::TwoPi()/2500);
  f1->SetParameter(6, 0);
  h2->Fit(f1);
  f1->SetNpx(100000);


  auto h3 = new TH1D("h3", "h3", 100, -1, 1);

  auto h4 = new TH2D("h4", "h4", numTimeBins, minTime, maxTime+1, numTimeBins, -30, 30);

  for(Long64_t entry=0; entry<tPat->GetEntries(); entry++){
    tPat->GetEntry(entry);

    if(pat->heading < -100) continue;

    Double_t thisHeading = pat->heading;
    if(thisHeading - lastHeading >= 180){
      while(thisHeading - lastHeading >= 180){
	thisHeading -= 360;
      }
    }
    if(thisHeading - lastHeading < -180){
      while(thisHeading - lastHeading < -180){
	thisHeading += 360;
      }
    }
    

    Double_t deltaHeading = thisHeading - lastHeading;    

    if(deltaHeading > 360){
      deltaHeading -= 360;
    }
    else if(deltaHeading < -360){
      deltaHeading += 360;
    }

    Double_t deltaHeadingExpected = f1->Eval(pat->realTime);
    
    if(TMath::Abs(deltaHeadingExpected - deltaHeading) < 0.1){
      h4->Fill(pat->realTime, deltaHeading);
    }
    h3->Fill(deltaHeadingExpected - deltaHeading);
    
    hRawHeading->Fill(pat->realTime, thisHeading);
    
    lastHeading = thisHeading;
  }

  auto c2 = new TCanvas();
  hRawHeading->Draw("colz");
  auto hRawHeading_px = hRawHeading->ProfileX();

  auto c3 = new TCanvas();
  h3->Draw("e");
  c3->SetLogy(1);

  auto c4 = new TCanvas();
  h4->Draw("colz");
  auto h5 = h4->ProfileX();
  h5->Draw("same");

  TF1* f2 = (TF1*) f1->Clone("f2");
  h5->Fit(f2);

  TGraph* gr = new TGraph();
  for(int binx=1; binx<=h5->GetNbinsX(); binx++){
    gr->SetPoint(gr->GetN(), h5->GetBinLowEdge(binx), h5->GetBinContent(binx));
  }
  gr->Draw();

  auto c5 = new TCanvas();
  auto gr2 = FancyFFTs::getPowerSpectrumTGraph(gr->GetN(), gr->GetY(), 1, PowSpecNorm::kSum, false);
  FFTWComplex* theFFT = FFTtools::doFFT(gr->GetN(), gr->GetY());
  gr2->Draw();

  for(int i=0; i < gr2->GetN(); i++){
    double powVal = gr2->GetY()[i];
    if(powVal <= 0.5){
      theFFT[i].re = 0;
      theFFT[i].im = 0;
    }
  }

  double* filteredGrVals = FFTtools::doInvFFT(gr->GetN(), theFFT);
  auto gr3 = new TGraph(gr->GetN(), gr->GetX(), filteredGrVals);
  new TCanvas();
  h5->Draw();
  gr3->SetLineColor(kBlue);
  gr3->Draw("lsame");

  TGraph* grHeading = new TGraph();
  grHeading->SetPoint(0, minTime, firstHeading);
  TGraph* grHeading2 = new TGraph();
  TGraph* grHeading3 = new TGraph();  
  UInt_t evalTime = minTime;

  const Double_t fitTimeMins = 5;
  std::vector<TF1*> fits;
  std::vector<TF1*> theFitInts;  
  Double_t startFitTime = minTime;
  Double_t endFitTime = minTime + fitTimeMins*60;
  Int_t numFits = 0;
  while(startFitTime < maxTime){

    fits.push_back(0);
    theFitInts.push_back(0);    
    auto theFit  = (TF1*) f2->Clone();
    fits.at(numFits) = theFit;
    theFit->SetRange(startFitTime, endFitTime);
    h5->Fit(theFit, "QR0");//, startFitTime, endFitTime);

    startFitTime += 0.5*fitTimeMins*60;
    endFitTime += 0.5*fitTimeMins*60;
    cout << theFit->GetChisquare()/theFit->GetNDF() << endl;

    new TCanvas();
    h5->Draw();
    theFit->Draw("lsame");
    new TCanvas();    
    TGraph *g = (TGraph*) theFit->DrawIntegral("al");
    auto theFitInt = (TF1*) f1Int->Clone();
    theFitInts.at(numFits) = theFitInt;
    for(int param=0; param<=6; param++){
      theFitInt->FixParameter(param, theFit->GetParameter(param));
      cout << param << "\t" << theFit->GetParameter(param) << "\t" << theFitInt->GetParameter(param) << endl;
    }
    Double_t thisMinTime, thisMaxTime;
    theFit->GetRange(thisMinTime, thisMaxTime);
    theFitInt->SetRange(thisMinTime, thisMaxTime);
    Double_t interceptMaybe = theFitInt->Eval(thisMinTime);
    cout << interceptMaybe << endl;

    
    theFitInt->SetParameter(7, -interceptMaybe - 10);
    hRawHeading_px->Fit(theFitInt, "RQ0");//, thisMinTime, thisMaxTime);
    
    theFitInt->Draw("l");
    hRawHeading->Draw("samecolz");    
    hRawHeading_px->Draw("same");
    // g->Draw("lsame");

    numFits++;

    while(evalTime < endFitTime - 0.33*fitTimeMins*60){
      Double_t lastHeading = grHeading->GetY()[grHeading->GetN()-1];
      Double_t deltaHeading = theFit->Eval(evalTime);
      grHeading->SetPoint(grHeading->GetN(), evalTime, lastHeading + deltaHeading);

      Double_t betterHeading = theFitInt->Eval(evalTime);
      while(betterHeading >= 360){
	betterHeading -=360;
      }
      while(betterHeading < 0){
	betterHeading +=360;
      }
      grHeading2->SetPoint(grHeading2->GetN(), evalTime, betterHeading);

      evalTime++;      
    }
    cout << UInt_t(startFitTime) << "\t" << evalTime << "\t" << UInt_t(endFitTime) << "\t" << numFits << endl;
  }
  // return;
  
  for(Long64_t entry=0; entry<tPat->GetEntries(); entry++){
    tPat->GetEntry(entry);

    if(pat->heading < -100) continue;

    Double_t betterHeading = grHeading2->Eval(pat->realTime);
    Double_t deltaBetterHeading = betterHeading - pat->heading;
    if(deltaBetterHeading >= 180){
      deltaBetterHeading -= 360;
    }
    if(deltaBetterHeading < -180){
      deltaBetterHeading += 360;
    }
    
    grHeading3->SetPoint(grHeading3->GetN(), pat->realTime, deltaBetterHeading);      
  }

  for(int i=0 ; i < grHeading->GetN(); i++){
    Double_t heading = grHeading->GetY()[i];
    while(heading >= 360){
      heading -= 360;
    }
    grHeading->GetY()[i] = heading;
  }
  new TCanvas();

  tPat->Draw("heading:realTime", "heading > -50", "colz");
  // grHeading->SetLineColor(kMagenta);
  grHeading->SetLineWidth(2);
  grHeading->Draw("lsame");
  grHeading2->Draw("lsame");  
  // return 0;
  
  TH2D* h0 = new TH2D("h0", "h0", 128, -10, 10, 128, -20, 20);
  TH2D* h01 = new TH2D("h01", "h01", numTimeBins, minTime, maxTime, 128, -20, 20);
  TH2D* h02 = new TH2D("h02", "h02", numTimeBins, minTime, maxTime, 128, -20, 20);
  
  // TH2D* h = new TH2D("h", "h", numTimeBins, -10, 10, 2, 0, 2);

  for(Long64_t entry=0; entry<t->GetEntries(); entry++){
    t->GetEntry(entry);
    tHeader->GetEntryWithIndex(eventNumber);
    tPat->GetEntryWithIndex(header->realTime);

    Double_t newHeading = grHeading->Eval(header->realTime);
    
    Double_t deltaHeading =  newHeading - pat->heading;
    h0->Fill(deltaPhiDeg, deltaHeading);
    h01->Fill(header->realTime, deltaHeading);
    h02->Fill(header->realTime, deltaPhiDeg);
  }
  auto c0c = new TCanvas();
  h0->Draw("colz");
  c0c->SetLogz(1);
  new TCanvas();
  h01->Draw("colz");
  auto h02_px = h02->ProfileX();
  h02_px->Draw("samehist");
  h02_px->SetLineColor(kRed);  
  new TCanvas();
  h02->Draw("colz");
  auto h01_px = h01->ProfileX();
  h01_px->SetLineColor(kMagenta);
  h01_px->SetLineWidth(2);  
  h01_px->Draw("samehist");

  TGraph* grDeltaHeadingDer = new TGraph();
  TGraph* grDeltaHeadingDer2 = new TGraph();  
  for(int binx=1; binx<h01_px->GetNbinsX(); binx++){
    Double_t der = h01_px->GetBinContent(binx+1) - h01_px->GetBinContent(binx);
    grDeltaHeadingDer->SetPoint(grDeltaHeadingDer->GetN(), h01_px->GetBinLowEdge(binx+1), der);
    grDeltaHeadingDer2->SetPoint(grDeltaHeadingDer2->GetN(), h01_px->GetBinLowEdge(binx+1), TMath::Abs(der));    
  }
  grDeltaHeadingDer->Draw("lsame");
  h01_px->Draw("samehist");
  
  new TCanvas();
  grDeltaHeadingDer2->Draw("alp");


  TH2D* h03 = new TH2D("h03", "h03", numTimeBins, minTime, maxTime+1, 128, -20, 20);
  TH2D* h04 = new TH2D("h04", "h04", numTimeBins, minTime, maxTime, 128, -20, 20);
  TH2D* h05 = new TH2D("h05", "h05", numTimeBins, minTime, maxTime, 128, -20, 20);
  TH2D* h06 = new TH2D("h06", "h06", numTimeBins, minTime, maxTime, 128, -20, 20);
  
  const Double_t deltaBadTime = 60;
  UInt_t lastBadTime = -1;
  for(Long64_t entry=0; entry<t->GetEntries(); entry++){
    t->GetEntry(entry);
    tHeader->GetEntryWithIndex(eventNumber);
    tPat->GetEntryWithIndex(header->realTime);

    if(pat->heading < -100) continue;

    if(grDeltaHeadingDer2->Eval(header->realTime) < 0.3){
      if(header->realTime - lastBadTime > deltaBadTime){
	h03->Fill(header->realTime, deltaPhiDeg);
      }
    }
    else{
      lastBadTime = header->realTime;
    }

    Double_t betterHeading = grHeading2->Eval(header->realTime);
    Double_t deltaBetterHeading = betterHeading - pat->heading;
    if(deltaBetterHeading >= 180){
      deltaBetterHeading -= 360;
    }
    if(deltaBetterHeading < -180){
      deltaBetterHeading += 360;
    }

    h04->Fill(header->realTime, deltaBetterHeading);
    h05->Fill(header->realTime, deltaPhiDeg+deltaBetterHeading);

    if(TMath::Abs(deltaBetterHeading) >= 2.5){

    // Int_t theBin = h5->FindBin(header->realTime);
    // // cout << theBin << "\t" << h5->GetBinError(theBin) << "\t" << h5->GetBinContent(theBin) << std::endl;
    // if(h5->GetBinError(theBin) > 0.3){
      // h06->Fill(header->realTime, deltaPhiDeg + deltaBetterHeading);
    }
    else{
      h06->Fill(header->realTime, deltaPhiDeg);
    }
  }

  auto c03 = new TCanvas();
  c03->Divide(2);
  c03->cd(1);
  h03->Draw("colz");
  c03->cd(2);
  h0->ProjectionX()->DrawNormalized("");      
  auto h03_px = h03->ProjectionY();
  h03_px->SetLineColor(kRed);
  h03_px->DrawNormalized("same");

  new TCanvas();
  h04->Draw("colz");
  grHeading3->Draw("lsame");

  new TCanvas();
  h05->Draw("colz");

  new TCanvas();
  h06->Draw("colz");

  new TCanvas();
  h02->Draw("colz");

  new TCanvas();
  auto h02_py = h02->ProjectionY();
  h02_py->SetLineColor(kRed);
  h02_py->DrawNormalized();
  h06->ProjectionY()->DrawNormalized("same");
  
}
