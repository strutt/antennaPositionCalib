{

  auto f = TFile::Open("fitPitchRollOffsetsPlots.root");

  TCanvas* c1 = new TCanvas();
  gPad->DrawFrame(0,-3,360,3);
  hDeltaPhiDeg2->Draw("colz");
  // hDeltaPhiDeg2->GetYaxis()->SetRangeUser(-3, 3);
  hDeltaPhiDeg2->SetTitle("#delta#phi vs #phi_{zoom} No dynamic correction");

  hDeltaPhiDeg2_pfx->Draw("same");
  hDeltaPhiDeg2_pfx->SetLineWidth(3);
  // hDeltaPhiDeg2_pfx->SetLineColor(6);  
  

  TCanvas* c2 = new TCanvas();
  gPad->DrawFrame(0,-3,360,3);
  hDeltaThetaDeg2->Draw("colz");
  hDeltaThetaDeg2->GetYaxis()->SetRangeUser(-3, 3);
  hDeltaThetaDeg2->SetTitle("#delta#theta vs #phi_{zoom} No dynamic correction");
  
  hDeltaThetaDeg2_pfx->Draw("same");
  hDeltaThetaDeg2_pfx->SetLineWidth(3);
  // hDeltaThetaDeg2_pfx->SetLineColor(6);      


  TCanvas* c3 = new TCanvas();
  gPad->DrawFrame(0,-3,360,3);  
  hDeltaPhiDeg2_pfx->SetTitle("#delta#phi vs #phi_{zoom} Fit summary plots");  
  hDeltaPhiDeg2_pfx->Draw();
  grDeltaPhi->SetLineColor(2);
  grDeltaPhi->Draw("lsame");
  grDeltaPhiNew->SetLineColor(4);
  grDeltaPhiNew->Draw("lsame");
  TLegend* l3 = new TLegend(0.15, 0.7, 0.55, 0.85);
  l3->AddEntry(hDeltaPhiDeg2_pfx, "Old expected - measured", "l");
  l3->AddEntry(grDeltaPhi, "New expected - old expected", "l");
  l3->AddEntry(grDeltaPhiNew, "Measured - new expected", "l");
  l3->Draw();
  c3->Update();

  TCanvas* c4 = new TCanvas();
  gPad->DrawFrame(0,-3,360,3);  
  hDeltaThetaDeg2_pfx->SetTitle("#delta#theta vs #phi_{zoom} Fit summary plots");
  hDeltaThetaDeg2_pfx->Draw();
  grDeltaTheta->SetLineColor(2);
  grDeltaTheta->Draw("lsame");
  grDeltaThetaNew->SetLineColor(4);
  grDeltaThetaNew->Draw("lsame");
  TLegend* l4 = new TLegend(0.15, 0.7, 0.55, 0.85);
  l4->AddEntry(hDeltaThetaDeg2_pfx, "Old expected - measured", "l");
  l4->AddEntry(grDeltaTheta, "New expected - old expected", "l");
  l4->AddEntry(grDeltaThetaNew, "Measured - new expected", "l");
  l4->Draw();
  c4->Update();

  
  
}
