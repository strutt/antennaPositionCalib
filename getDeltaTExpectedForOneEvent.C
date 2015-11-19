{
  gSystem->Load("libAnitaCorrelator.so");

  TFile* f = TFile::Open("~/UCL/ANITA/flight1415/root/run352/gpsFile352.root");
  TTree* t = (TTree*) f->Get("adu5PatTree");

  Adu5Pat* pat = NULL;
  t->SetBranchAddress("pat", &pat);
  
  t->BuildIndex("realTime");
  t->GetEntryWithIndex(1420170384);
  
  UsefulAdu5Pat usefulPat(pat);

  cout << usefulPat.getDeltaTExpected(16, 0, 315*TMath::DegToRad(), 6*TMath::DegToRad()) << endl;





}
