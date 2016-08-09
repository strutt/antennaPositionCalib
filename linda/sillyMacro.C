void sillyMacro(){

  char headerName[1000];
  char gpsName[1000];
  char corName[1000];

  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  Adu5Pat *pat2 =0;
  CorrelationSummaryAnita3 *cor=0;

  TChain *gpsChain = new TChain("adu5PatTree");
  TChain *gpsChain2 = new TChain("adu5PatTree");
  TChain *headChain = new TChain("headTree");
  TChain *corChain = new TChain("corTree");

  for (unsigned int run=331;run<355;++run){
    sprintf(headerName,"/unix/anita3/flight1415/root/run%d/timedHeadFile%d.root",run,run);
    sprintf(gpsName,"/unix/anita3/flight1415/root/run%d/gpsEvent%d.root",run,run);
    sprintf(corName, "/unix/anita3/linda/corTrees/corRun_NEW11_HPOL_%d.root", run);
    headChain->Add(headerName);
    gpsChain->Add(gpsName);
    gpsChain2->Add(gpsName);
    corChain->Add(corName);
  }
  headChain->SetBranchAddress("header",&header);
  gpsChain->SetBranchAddress("pat",&pat);
  gpsChain2->SetBranchAddress("pat",&pat2);
  corChain->SetBranchAddress("cor",&cor);

  UInt_t gps_eventNumber;
  gpsChain->SetBranchAddress("eventNumber", &gps_eventNumber);
  UInt_t gps2_eventNumber;
  gpsChain2->SetBranchAddress("eventNumber", &gps2_eventNumber);

  headChain->BuildIndex("header->eventNumber");
  gpsChain->BuildIndex("eventNumber");
  gpsChain2->BuildIndex("pat->realTime");


  int maxEntry=corChain->GetEntries();

  // AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;                                                    
  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  char cpol[100];

  Double_t phiWave, thetaWave, lower, upper;
  Double_t maxCorrTime, deltaTExpected;
  int ant1, ant2;
  headChain->GetEntry(0);
  double firstTS = header->triggerTime;
  headChain->GetEntry(headChain->GetEntries()-1);
  double lastTS = header->triggerTime;

  Double_t additionalPhi = 22.5*TMath::DegToRad();
  Double_t TwoPi = TMath::Pi()*2.;
  //  TH2D *h1 = new TH2D("h1", "", 100, firstTS, lastTS, 100, 0, 3600*24);                               
  TH2D *h1 = new TH2D("h1", "", 100, firstTS, lastTS, 100, -49.5, 50.5);
  TH1D *dHeading = new TH1D("dHeading", "", 100, -5, 5);

  for(Long64_t entry=0;entry<maxEntry;++entry) {
    corChain->GetEntry(entry);

    Long64_t headEntry = headChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(headEntry < 0 ) continue;
    headChain->GetEntry(headEntry);

    Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(cor->eventNumber);
    if(gpsEntry < 0 ) continue;
    gpsChain->GetEntry(gpsEntry);

    Long64_t gpsEntry2 = gpsChain2->GetEntryNumberWithIndex(header->triggerTime);
    if(gpsEntry2 < 0 ) continue;
    gpsChain2->GetEntry(gpsEntry2);

    if(gps_eventNumber!=gps2_eventNumber) cout << header->eventNumber <<  " " << gps2_eventNumber << endl;

    //    if ((header->triggerTime%(3600*24))!=(pat->timeOfDay/1e3)) cout << header->triggerTime%(3600*24) << " " << pat->timeOfDay/1e3 << endl;
    double tday = (header->triggerTime%(3600*24));
    h1->Fill(header->triggerTime, (tday - pat->timeOfDay/1e3) );

    dHeading->Fill((pat->heading-pat2->heading));

  }

  h1->Draw("colz");


  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(1);
  dHeading->SetTitle("Heading(eventNumber)-heading(realTime);Heading(eventNumber)-heading(realTime) [degrees];Number of events");
  dHeading->Draw();
  c1->Print("dheading.png");

}
