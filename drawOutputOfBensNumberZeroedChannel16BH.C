void drawOutputOfBensNumberZeroedChannel16BH(){

  // gStyle->SetOptStat("mre");
    

  
  auto cOld = new TChain("angResTree");
  // cOld->Add("januaryTests/photogrammetryNumbers/generateAngularResolutionTreePlots_3*.root");
  cOld->Add("januaryTests/photogrammetryNumbersZeroedChannel16BH_cosminV3Trees/generateAngularResolutionTreePlots_3*.root");  

  auto cNew = new TChain("angResTree");
  cNew->Add("januaryTests/bensNumbersZeroedChannel16BH_cosminV3Trees/generateAngularResolutionTreePlots_3*.root");
  
  const Int_t numVars = 2;
  TString drawVars[numVars] = {"deltaPhiDeg", "deltaThetaDeg"};
  TString varNames[numVars] = {"#delta#phi", "#delta#theta"};  

  for(int varInd=0; varInd < numVars; varInd++){
  
    auto c1 = new TCanvas();
    const Int_t numPhiBins = 128;

    TProfile* hProfOld = new TProfile("hProfOld", "Photogrammetry", numPhiBins, 0, 360);
    TH1D* hOld = new TH1D("hOld", "Photogrammetry", numPhiBins/2, -5, 5);    

    TString drawCommandOld = TString::Format("%s:phiExpected>>hProfOld", drawVars[varInd].Data());
    cOld->Draw(drawCommandOld, "TMath::Abs(deltaPhiDeg) < 5", "goff");

    TString drawCommandOld2 = TString::Format("%s>>hOld", drawVars[varInd].Data());
    cOld->Draw(drawCommandOld2, "TMath::Abs(deltaPhiDeg) < 5", "goff");    

    hProfOld->SetLineColor(kBlue);
    hOld->SetLineColor(kBlue);    


    TProfile* hProfNew = new TProfile("hProfNew", "Ben's numbers", numPhiBins, 0, 360);
    TH1D* hNew = new TH1D("hNew", "Ben's numbers", numPhiBins/2, -5, 5);

    TString drawCommandNew = TString::Format("%s:phiExpected>>hProfNew", drawVars[varInd].Data());
    cNew->Draw(drawCommandNew, "TMath::Abs(deltaPhiDeg) < 5", "goff");
    TString drawCommandNew2 = TString::Format("%s>>hNew", drawVars[varInd].Data());    
    cNew->Draw(drawCommandNew2, "TMath::Abs(deltaPhiDeg) < 5", "goff");    

    hProfNew->SetLineColor(kRed);
    hNew->SetLineColor(kRed);    
  
    hProfNew->Draw();  
    hProfOld->Draw("same");

    // hProfOld->Draw();
    // hProfNew->Draw("same");  

    auto l = c1->BuildLegend();
    l->Draw();


    // TH1D* hs[2] = {hOld, hNew};
    TH1D* hs[2] = {hNew, hOld};    
    TCanvas* c2 = RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");
    auto l2 = c2->BuildLegend();
    l2->Draw();
    
    hProfNew->SetTitle("Comparison of WAIS divide photogrammetry vs. Ben's numbers "+varNames[varInd]+" ; Expected #phi (Degrees); " + varNames[varInd] + " (Degrees)");
    hProfOld->SetTitle("Comparison of WAIS divide photogrammetry vs. Ben's numbers "+varNames[varInd]+" ; Expected #phi (Degrees); " + varNames[varInd] + " (Degrees)");    
    c1->Update();

    hOld->SetTitle("Comparison of WAIS divide pulses "+varNames[varInd]+", photogrammetry vs. Ben's numbers;"+varNames[varInd]+" (Degrees); Events/bin");
    hNew->SetTitle("Comparison of WAIS divide pulses "+varNames[varInd]+", photogrammetry vs. Ben's numbers;"+varNames[varInd]+" (Degrees); Events/bin");    
    c2->Update();
  }
}
