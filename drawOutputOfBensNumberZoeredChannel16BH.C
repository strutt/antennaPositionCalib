
void drawOutputOfTestingZeroedChannel16BH(){

  gStyle->SetOptStat("mre");
    

  
  auto fOld = TFile::Open("januaryTests/photogrammetryNumbers/generateAngularResolutionTreePlots_352_2016-01-07_15-14-59.root");

  const Int_t numVars = 2;
  TString drawVars[numVars] = {"deltaPhiDeg", "deltaThetaDeg"};
  TString varNames[numVars] = {"#delta#phi", "#delta#theta"};  

  for(int varInd=0; varInd < numVars; varInd++){
  
    auto c1 = new TCanvas();
    const Int_t numPhiBins = 128;

    TProfile* hProfOld = new TProfile("hProfOld", "Including 16BH", numPhiBins, 0, 360);
    TH1D* hOld = new TH1D("hOld", "Including 16BH", numPhiBins/2, -5, 5);    
    auto tOld = (TTree*) fOld->Get("angResTree");

    TString drawCommandOld = TString::Format("%s:phiExpected>>hProfOld", drawVars[varInd].Data());
    tOld->Draw(drawCommandOld, "TMath::Abs(deltaPhiDeg) < 5", "goff");

    TString drawCommandOld2 = TString::Format("%s>>hOld", drawVars[varInd].Data());
    tOld->Draw(drawCommandOld2, "TMath::Abs(deltaPhiDeg) < 5", "goff");    

    hProfOld->SetLineColor(kBlue);
    hOld->SetLineColor(kBlue);    

    auto fNew = TFile::Open("generateAngularResolutionTreePlots_352_2016-01-29_12-49-06.root");
    TProfile* hProfNew = new TProfile("hProfNew", "Zeroing 16BH", numPhiBins, 0, 360);
    TH1D* hNew = new TH1D("hNew", "Zeroing 16BH", numPhiBins/2, -5, 5);
    auto tNew = (TTree*) fNew->Get("angResTree");

    TString drawCommandNew = TString::Format("%s:phiExpected>>hProfNew", drawVars[varInd].Data());
    tNew->Draw(drawCommandNew, "TMath::Abs(deltaPhiDeg) < 5", "goff");
    TString drawCommandNew2 = TString::Format("%s>>hNew", drawVars[varInd].Data());    
    tNew->Draw(drawCommandNew2, "TMath::Abs(deltaPhiDeg) < 5", "goff");    

    hProfNew->SetLineColor(kRed);
    hNew->SetLineColor(kRed);    
  

    hProfOld->Draw();
    hProfNew->Draw("same");  

    auto l = c1->BuildLegend();
    l->Draw();


    TH1D* hs[2] = {hNew, hOld};
    TCanvas* c2 = RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");
    auto l2 = c2->BuildLegend();
    l2->Draw();
    
    hProfOld->SetTitle("Comparison of WAIS divide photogrammetry "+varNames[varInd]+" when channel 16BH is zeroed (run 352 only); Expected #phi (Degrees); " + varNames[varInd] + " (Degrees)");
    c1->Update();

    hNew->SetTitle("Comparison of WAIS divide photogrammetry "+varNames[varInd]+" when channel 16BH is zeroed (run 352 only); " + varNames[varInd] + " (Degrees); Events/bin");
    c2->Update();
  }
}
