void drawMyVPolFitterImplementationPlots(){

  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-17-57.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-51-19.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_12-08-13.root");

  TString fName = "myVPolFitterImplementationPlots_2016-04-13_14-28-17.root";
  TFile* f = TFile::Open(fName);

  EColor cols[NUM_PHI] = {kRed, kBlue, kGreen, kOrange, kCyan, kMagenta, kYellow, kBlack,
			  kRed, kBlue, kGreen, kOrange, kCyan, kMagenta, kYellow, kBlack}; 


  const int numFitSteps = 2;
  for(int fitInd=0; fitInd < numFitSteps; fitInd++){
    TString canName = TString::Format("iteration_%d", fitInd);
    auto c = new TCanvas(canName, canName, 1600, 1600);
    c->Divide(2, 3);
    for(int vertAdj=0; vertAdj<2; vertAdj++){
      TString keyWord = vertAdj == 0 ? "Adjacent" : "Vertical";
      for(int ring=0; ring < AnitaRing::kNotARing; ring++){
	c->cd(ring*2 + vertAdj + 1);
	for(int phi=0; phi < NUM_PHI; phi++){
	  int ant = NUM_RING*ring + phi;
	  TString name = TString::Format("gr%sDists_%d_%d", keyWord.Data(), ant, fitInd);
	  TGraph* gr = (TGraph*) f->Get(name);
	  gr->SetMarkerColor(cols[phi]);
	  gr->SetMarkerStyle(1);	
	  TString drawOpt = phi==0 ? "ap" : "psame";
	  gr->Draw(drawOpt);
	  if(phi==0){
	    TString title = "#deltat_{measured} - #deltat_{expected} " + keyWord + " pairs ";
	    switch(ring){
	    case 0:
	      title += "top ring ";
	      break;
	    case 1:
	      title += "middle ring ";
	      break;
	    case 2:
	      title += "bottom ring ";
	      break;
	    }
	    title += TString::Format(" fit iteration %d", fitInd);
	    title += "; #phi_{wave} (Radians); #Delta#deltat (ns)";
	    gr->SetTitle(title);
	    gr->GetXaxis()->SetLimits(0, TMath::TwoPi());
	    gr->GetYaxis()->SetRangeUser(-1, 1);
	  }
	}
      }
    }
  }
  // return;
  
  TGraph* grAdjMean_0 = (TGraph*) f->Get("grAdjMean_0");
  TGraph* grVertMean_0 = (TGraph*) f->Get("grVertMean_0");
  TGraph* grAdjRms_0 = (TGraph*) f->Get("grAdjRms_0");
  TGraph* grVertRms_0 = (TGraph*) f->Get("grVertRms_0");
  TGraph* grAdjSlope_0 = (TGraph*) f->Get("grAdjSlope_0");
  TGraph* grVertSlope_0 = (TGraph*) f->Get("grVertSlope_0");
  TGraph* grAdjMean_1 = (TGraph*) f->Get("grAdjMean_1");
  TGraph* grVertMean_1 = (TGraph*) f->Get("grVertMean_1");
  TGraph* grAdjRms_1 = (TGraph*) f->Get("grAdjRms_1");
  TGraph* grVertRms_1 = (TGraph*) f->Get("grVertRms_1");
  TGraph* grAdjSlope_1 = (TGraph*) f->Get("grAdjSlope_1");
  TGraph* grVertSlope_1 = (TGraph*) f->Get("grVertSlope_1");


  TGraph* grAdjMean_0_weighted = (TGraph*) f->Get("grAdjMean_0_weighted");
  TGraph* grVertMean_0_weighted = (TGraph*) f->Get("grVertMean_0_weighted");
  TGraph* grAdjRms_0_weighted = (TGraph*) f->Get("grAdjRms_0_weighted");
  TGraph* grVertRms_0_weighted = (TGraph*) f->Get("grVertRms_0_weighted");
  TGraph* grAdjSlope_0_weighted = (TGraph*) f->Get("grAdjSlope_0_weighted");
  TGraph* grVertSlope_0_weighted = (TGraph*) f->Get("grVertSlope_0_weighted");
  TGraph* grAdjMean_1_weighted = (TGraph*) f->Get("grAdjMean_1_weighted");
  TGraph* grVertMean_1_weighted = (TGraph*) f->Get("grVertMean_1_weighted");
  TGraph* grAdjRms_1_weighted = (TGraph*) f->Get("grAdjRms_1_weighted");
  TGraph* grVertRms_1_weighted = (TGraph*) f->Get("grVertRms_1_weighted");
  TGraph* grAdjSlope_1_weighted = (TGraph*) f->Get("grAdjSlope_1_weighted");
  TGraph* grVertSlope_1_weighted = (TGraph*) f->Get("grVertSlope_1_weighted");
  TGraph* grTotal_0 = (TGraph*) f->Get("grTotal_0");
  TGraph* grTotal_1 = (TGraph*) f->Get("grTotal_1");
    

  auto c0 = new TCanvas();
  grTotal_0->Draw();
  grTotal_0->SetLineColor(kBlack);  
  grAdjMean_0_weighted->Draw("lsame");
  grAdjMean_0_weighted->SetLineColor(kBlue);
  grVertMean_0_weighted->Draw("lsame");
  grVertMean_0_weighted->SetLineColor(kRed);
  grAdjRms_0_weighted->Draw("lsame");
  grAdjRms_0_weighted->SetLineColor(kMagenta);
  grVertRms_0_weighted->Draw("lsame");
  grVertRms_0_weighted->SetLineColor(kOrange);
  grAdjSlope_0_weighted->Draw("lsame");
  grAdjSlope_0_weighted->SetLineColor(kCyan);
  grVertSlope_0_weighted->Draw("lsame");
  grVertSlope_0_weighted->SetLineColor(kSpring);
  c0->BuildLegend();

  auto c01 = new TCanvas();
  grTotal_1->Draw();
  grTotal_1->SetLineColor(kBlack);  
  grAdjMean_1_weighted->Draw("lsame");
  grAdjMean_1_weighted->SetLineColor(kBlue);
  grVertMean_1_weighted->Draw("lsame");
  grVertMean_1_weighted->SetLineColor(kRed);
  grAdjRms_1_weighted->Draw("lsame");
  grAdjRms_1_weighted->SetLineColor(kMagenta);
  grVertRms_1_weighted->Draw("lsame");
  grVertRms_1_weighted->SetLineColor(kOrange);
  grAdjSlope_1_weighted->Draw("lsame");
  grAdjSlope_1_weighted->SetLineColor(kCyan);
  grVertSlope_1_weighted->Draw("lsame");
  grVertSlope_1_weighted->SetLineColor(kSpring);
  c01->BuildLegend();

  
  auto c1 = new TCanvas();
  grAdjMean_0->Draw("al");
  grAdjMean_1->Draw("lsame");
  grAdjMean_0->SetLineColor(kBlue);
  grAdjMean_1->SetLineColor(kRed);
  c1->BuildLegend();

  auto c2 = new TCanvas();
  grVertMean_0->Draw("al");
  grVertMean_1->Draw("lsame");
  grVertMean_0->SetLineColor(kBlue);
  grVertMean_1->SetLineColor(kRed);
  c2->BuildLegend();

  auto c3 = new TCanvas();
  grAdjRms_0->Draw("al");
  grAdjRms_1->Draw("lsame");
  grAdjRms_0->SetLineColor(kBlue);
  grAdjRms_1->SetLineColor(kRed);
  c3->BuildLegend();
  
  auto c4 = new TCanvas();
  grVertRms_0->Draw("al");
  grVertRms_1->Draw("lsame");
  grVertRms_0->SetLineColor(kBlue);
  grVertRms_1->SetLineColor(kRed);
  c4->BuildLegend();
  
  auto c5 = new TCanvas();
  grAdjSlope_0->Draw("al");
  grAdjSlope_1->Draw("lsame");
  grAdjSlope_0->SetLineColor(kBlue);
  grAdjSlope_1->SetLineColor(kRed);
  c5->BuildLegend();
  
  auto c6 = new TCanvas();
  grVertSlope_0->Draw("al");
  grVertSlope_1->Draw("lsame");
  grVertSlope_0->SetLineColor(kBlue);
  grVertSlope_1->SetLineColor(kRed);
  c6->BuildLegend();
  
  
  
  
}
