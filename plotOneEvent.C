std::vector<TGraph*> getDirectionTGraphs(Double_t heading, Int_t numDirs, TLegend* l){

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  Double_t northDeg = heading;
  const Double_t THETA_RANGE = 150;

  std::vector<TGraph*> grs;
  
  Double_t deltaPhiDeg = 360./numDirs;
  for(int dir=0; dir<numDirs; dir++){
    Double_t ys[2] = {THETA_RANGE/2, -THETA_RANGE/2};
    Double_t phiDeg = northDeg - dir*deltaPhiDeg;
    if(phiDeg < 0) phiDeg += 360;
    Double_t xs[2] = {phiDeg, phiDeg};
    grs.push_back(new TGraph(2, xs, ys));
    // grs.at(dir)->SetLineColor(kWhite);
    grs.at(dir)->SetLineColor(kGray);    
    grs.at(dir)->SetLineStyle(1+dir);
    grs.at(dir)->SetLineWidth(2);

  }

  l->AddEntry(grs.at(0), "North", "l");
  l->AddEntry(grs.at(1), "West", "l");
  l->AddEntry(grs.at(2), "South", "l");
  l->AddEntry(grs.at(3), "East", "l");  
  //  l->SetFillColor(kGray);

  return grs;

}





void plotOneEvent(){

  // Load libraries into ROOT
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->Load("libBensAnitaTools.so");  
  gSystem->Load("libAnitaEvent.so");
  
  TString calEventFilePath = "~/UCL/ANITA/flight1415/root/run352/calEventFile352.root";
  TString headFilePath = "~/UCL/ANITA/flight1415/root/run352/headFile352.root";
  
  TFile* fE = TFile::Open(calEventFilePath);
  TTree* tE = (TTree*) fE->Get("eventTree");

  CalibratedAnitaEvent* calEvent = NULL;
  tE->SetBranchAddress("event", &calEvent);  
  // RawAnitaEvent* rawEvent = NULL;
  // tE->SetBranchAddress("event", &rawEvent);  

  TFile* fH = TFile::Open(headFilePath);
  TTree* tH = (TTree*) fH->Get("headTree");
  RawAnitaHeader* header = NULL;
  tH->SetBranchAddress("header", &header);

  Long64_t numEntries = tE->GetEntries();


  const Long64_t entryIWant = 496;
  tH->GetEntry(entryIWant);
  tE->GetEntry(entryIWant);

  // UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(rawEvent, WaveCalType::kDefault, header);
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(calEvent);//, WaveCalType::kDefault, header);

  // CrossCorrelator* cc = NULL;
  CrossCorrelator* cc = new CrossCorrelator();  
  cc->correlateEvent(realEvent, AnitaPol::kHorizontal);
  
  Double_t a, b, c;
  TH2D* hImage = cc->makeTriggeredImage(AnitaPol::kHorizontal, a, b, c, header->l3TrigPatternH);
  hImage->Draw("colz");

  Double_t phiMin = -11.25-45;//cc->getBin0PhiDeg();
  Double_t phiMax = phiMin + 360;
  Double_t thetaMin = -150/2;
  Double_t thetaMax = 150/2;

  TH2D* hBlank = new TH2D("hBlank", "asd",
  			  16*64, phiMin, phiMax,
  			  128, thetaMin, thetaMax);
  hBlank->GetXaxis()->SetTitle("Azimuth (degrees)");
  hBlank->GetXaxis()->SetTitle("");  
  hBlank->GetXaxis()->SetNdivisions(0);  

  // for(int phi=0; phi<16; phi++){
  //   hBlank->GetXaxis()->SetBinLabel(phi+1, TString::Format("%d", phi+1));
  // }
  // hBlank->GetXaxis()->SetLabelOffset(0.005);
  hBlank->GetYaxis()->SetTitle("Elevation (Degrees)");
  hBlank->Draw();
  // hBlank->GetXaxis()->SetNdivisions(16);
  
  // hImage->Draw("colzsame");
  hBlank->GetYaxis()->Draw();
  TF1 *f1=new TF1("fAxisDrawFunction","x",0.5,16.5);
  TGaxis *A1 = new TGaxis(phiMin,thetaMin,phiMax,thetaMin,"fAxisDrawFunction",16,"+");
  A1->SetTitle("Phi-sector");
  A1->SetLabelSize(0.03);
  A1->SetTickSize(0);
  A1->SetTickLength(0);
  
  A1->Draw();
  
  delete f1;
  gPad->RedrawAxis();

  Double_t northDeg = 300;
  const Int_t numDirs = 4;
  TLegend* l0 = new TLegend(0.8, 0.8, 1, 1);
  std::vector<TGraph*> grs = getDirectionTGraphs(northDeg, numDirs, l0);
  for(UInt_t i=0; i<grs.size(); i++){
    grs[i]->Draw("lsame");
  }
  l0->Draw();


  return;
  
  new TCanvas();  
  TH2D* hImage2 = cc->makeZoomedImage(AnitaPol::kHorizontal, header->l3TrigPatternH, b, c);  
  hImage2->Draw("colz");
  


  // new TCanvas();

  // TGraph* gr = realEvent->getGraph(0, AnitaPol::kHorizontal);
  // TGraph* gr2 = realEvent->getGraph(16, AnitaPol::kHorizontal);  
  // TGraph* gr3 = cc->getCrossCorrelationGraph(AnitaPol::kHorizontal, 0, 16);

  
  // gr->Draw();
  // gr2->SetLineColor(kRed);
  // gr2->Draw("lsame");
  // new TCanvas();
  // gr3->Draw();

  

  delete realEvent;


}
