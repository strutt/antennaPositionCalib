{

  TFile* f = TFile::Open("compareMeasuredWithExpectedDeltaTsPlots.root");


   const int NRGBs = 3, NCont = 999;
   gStyle->SetNumberContours(NCont);
   Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
   Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
   Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
   Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont); 

  
  // // Draw feed and photogrammetry examples for 0, 16.
  // TH2D* hPureDiff_2D_feed_0_16 = (TH2D*) f->Get("hPureDiff_2D_feed_0_16");
  // TH2D* hPureDiff_2D_photo_0_16 = (TH2D*) f->Get("hPureDiff_2D_photo_0_16");

  // TCanvas* c1 = new TCanvas();
  // hPureDiff_2D_feed_0_16->Draw("colz");
  // hPureDiff_2D_feed_0_16->GetXaxis()->SetRangeUser(250, 360);
  // hPureDiff_2D_feed_0_16->GetYaxis()->SetRangeUser(-8, -4);

  // hPureDiff_2D_feed_0_16->SetMaximum(0.25);
  // hPureDiff_2D_feed_0_16->SetMinimum(-0.25);

  // TCanvas* c2 = new TCanvas();
  // hPureDiff_2D_photo_0_16->Draw("colz");
  // hPureDiff_2D_photo_0_16->GetXaxis()->SetRangeUser(250, 360);
  // hPureDiff_2D_photo_0_16->GetYaxis()->SetRangeUser(-8, -4);

  // hPureDiff_2D_photo_0_16->SetMaximum(0.25);
  // hPureDiff_2D_photo_0_16->SetMinimum(-0.25);


  // TCanvas* c3 = new TCanvas();
  // TH1D* hPureDiff_feed_0_16 = (TH1D*) f->Get("hPureDiff_feed_0_16");
  // TH1D* hPureDiff_photo_0_16 = (TH1D*) f->Get("hPureDiff_photo_0_16");
  // hPureDiff_feed_0_16->Draw();
  // hPureDiff_feed_0_16->SetLineColor(kRed);
  // hPureDiff_photo_0_16->Draw("same");
  // hPureDiff_photo_0_16->SetLineColor(kBlue);

  // hPureDiff_feed_0_16->SetTitle("Distributions of #deltat_{expected} - #deltat_{measured}");
  // TH1D* hs[2] = {hPureDiff_feed_0_16, hPureDiff_photo_0_16};
  // TString titles[2] = {"Feed locations", "Photogrammetry locations"};
  // TLegend* l3 = RootTools::makeLegend(hs, 2, titles, "l", 0.6, 0.65, 0.95, 0.85);
  // l3->Draw();
  


  TCanvas* c4 = new TCanvas();
  TH1D* hPureDiff_feed_0_16 = (TH1D*) f->Get("hPureDiff_feed_0_16");
  TH1D* hPureDiff_photo_0_16 = (TH1D*) f->Get("hPureDiff_photo_0_16");
  hPureDiff_feed_0_16->Draw();
  hPureDiff_feed_0_16->SetLineColor(kRed);
  hPureDiff_photo_0_16->Draw("same");
  hPureDiff_photo_0_16->SetLineColor(kBlue);

  hPureDiff_feed_0_16->SetTitle("Distributions of #deltat_{expected} - #deltat_{measured}");
  TH1D* hs[2] = {hPureDiff_feed_0_16, hPureDiff_photo_0_16};
  TString titles[2] = {"Feed locations", "Photogrammetry locations"};
  TLegend* l4 = RootTools::makeLegend(hs, 2, titles, "l", 0.6, 0.65, 0.95, 0.85);
  l4->Draw();

  
  // RootTools::makeZaxisScaleEqualAboutZero(hPureDiff_2D_feed_0_16);
  TCanvas* c5 = new TCanvas();
  TGraph* grDeltaRScan= (TGraph*) f->Get("grDeltaRScan");
  TGraph* grFitterMin = (TGraph*) f->Get("grFitterMin");
  
  grDeltaRScan->Draw("al");
  grFitterMin->Draw("psame");  


  // grDeltaRScan->SetMinimum(0);
  grDeltaRScan->GetXaxis()->SetNoExponent(1);  
  
  grDeltaRScan->SetLineColor(kBlack);
  grDeltaRScan->SetMarkerStyle(0);  
  grFitterMin->SetMarkerStyle(8);
  grFitterMin->SetMarkerColor(kRed);
  TLegend* l5 = new TLegend(0.65, 0.8, 1, 1);
  l5->AddEntry(grDeltaRScan, "Manual scan in #deltar", "l");
  l5->AddEntry(grFitterMin, "Minuit2 minimum", "p");
  l5->Draw();

  
  return;

}
