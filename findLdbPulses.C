
    // \item LDB Discone (10\,kV) (runs 145-149), delay 25, 225, 425, 625, 825 ms;
    // \item LDB Seavey (6\,kV) (runs 151-152), delay 25, 225, 425, 625, 825 ms;
    // \item LDB Seavey (10\,kV) (runs 154-161), delay 50, 250, 450, 650, 850 ms;
{

  auto c = new TChain("headTree");
  for(Int_t run=145; run <= 161; run++){
    auto fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    c->Add(fileName);
  }
  RawAnitaHeader* header= NULL;
  c->SetBranchAddress("header", &header);
  c->Show(0);
  c->Show(c->GetEntries()-1);
  Double_t repeatTime = 2e8;
  auto h2 = new TH2D("h2", "VPOL selection; triggerTimeNs (ns); eventNumber",  200, 0, repeatTime, /*1076*/ 2657, 8165401, 10823328); //9241486);
  TCanvas* c1 = new TCanvas();
  c->Draw("eventNumber:(triggerTimeNs-25)%" + TString::Format("%lf>>h2", repeatTime),
	  "(trigType & 1)==1", "colz");


  TCanvas* c2 = new TCanvas();
  auto h = new TH1D("h", "VPOL selection; triggerTimeNs (ns)",  200, 0, repeatTime); //, /*1076*/ 2657, 8165401, 10823328); //9241486);
  c->Draw("(triggerTimeNs-25)%" + TString::Format("%lf>>h", repeatTime),
	  "(trigType & 1)==1", "");

  h->GetYaxis()->SetNoExponent(1);
  h->GetXaxis()->SetNoExponent(1);  


  h2->GetYaxis()->SetNoExponent(1);
  h2->GetXaxis()->SetNoExponent(1);  
  c1->SetLogz(1);
  


}
