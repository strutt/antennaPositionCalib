double findSlope(Int_t wrappedPhi, TGraph* gr){

  Double_t SUMx = 0;
  Double_t SUMy = 0;
  Double_t SUMxy = 0;
  Double_t SUMxx = 0;
   
  const int n = gr->GetN();
  for(int i=0; i < n; ++i){
    double thisX = gr->GetX()[i];
      
    if(wrappedPhi){
      if(thisX > TMath::Pi()){
	thisX -= TMath::TwoPi();
      }
    }

    SUMx = SUMx + thisX; //xIn[i];
    SUMy = SUMy + gr->GetY()[i];
    SUMxy = SUMxy + thisX*gr->GetY()[i];
    SUMxx = SUMxx + thisX*thisX;

  }
  Double_t slope = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx );

  return slope;  
}


TGraph* makeTestCase(){

  const int n = 100;
  Double_t xMax = 10;
  Double_t norm = xMax*xMax;
  Double_t dx = 2*xMax/n;
  TRandom3 rndm;
  TGraph* gr = new TGraph();
  for(int i=0; i<=n; i++){
    Double_t x = -xMax + dx*i;
    Double_t y = (10*x*x + rndm.Gaus(0, 1))/norm;
    gr->SetPoint(i, x, y);
  }

  gr->SetTitle("Simple quadratic curve plus small gaussian noise; x (no units); y (no units)");
  gr->Fit("pol1");
  gr->SetMarkerStyle(2);
  
  return gr;
}

void drawMyVPolFitterImplementationPlots(){

  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-17-57.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-51-19.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_12-08-13.root");

  // TString fName = "myVPolFitterImplementationPlots_2016-04-13_14-28-17.root";
  //  TString fName = "myVPolFitterImplementationPlots_2016-04-13_15-10-19.root"; // m10000r100s10000
  // TString fName = "myVPolFitterImplementationPlots_2016-04-13_17-52-53.root"; // m10000r10000s10000


  //  TString fName = "myVPolFitterImplementationPlots_2016-04-14_00-52-36.root"; // first t, (phi?), combined fom.
  
  // TString fName = "myVPolFitterImplementationPlots_2016-04-15_10-48-02.root"; // now really does Linda's fitting order

  // TString fName = "myVPolFitterImplementationPlots_2016-04-15_12-25-38.root";// phi (grad), maxCalls

  TString fName = "myVPolFitterImplementationPlots_2016-04-15_13-07-44.root "; //r, phi (grad, rms), maxCalls

  TFile* f = TFile::Open(fName);

  EColor cols[NUM_PHI] = {kRed, kBlue, kMagenta, kCyan, kGreen, kOrange, kViolet, kBlack,
			  kRed, kBlue, kMagenta, kCyan, kGreen, kOrange, kViolet, kBlack}; 

  
  if(true){
    const int numFitSteps = 2; // including iteration 0.
    for(int fitInd=0; fitInd < numFitSteps; fitInd++){
      TString canName = TString::Format("iteration_%d", fitInd);
      auto c = new TCanvas(canName, canName, 1600, 1600);
      c->Divide(2, 3);
      for(int vertAdj=0; vertAdj<2; vertAdj++){
	TString keyWord = vertAdj == 0 ? "Adjacent" : "Vertical";
	for(int ring=0; ring < AnitaRing::kNotARing; ring++){
	  c->cd(ring*2 + vertAdj + 1);
	  for(int phi=0; phi < NUM_PHI; phi++){
	    int ant = NUM_PHI*ring + phi;
	    TString name = TString::Format("gr%sDists_%d_%d", keyWord.Data(), ant, fitInd);
	    TGraph* gr = (TGraph*) f->Get(name);	 
	  
	    TProfile* p = new TProfile(TString::Format("p_%s", name.Data()), "", 128, 0, TMath::TwoPi());
	    // std::cout << fitInd << "\t" << vertAdj << "\t" << ring << "\t" << phi << "\t"
	    // 	    << p->GetName() << std::endl;
	    p->SetLineColor(cols[phi]);
	    gr->SetMarkerColor(cols[phi]);
	    gr->SetMarkerStyle(1);

	    for(int i=0; i< gr->GetN(); i++){
	      p->Fill(gr->GetX()[i], gr->GetY()[i]);
	    }
	    // TString drawOpt1 = phi==0 ? "ap" : "psame";
	    // gr->Draw(drawOpt1);

	    TString drawOpt2 = phi==0 ? "" : "lsame";	  
	    p->Draw(drawOpt2);	  
	    if(phi==0){
	      TString title = "#deltat_{measured} - #deltat_{expected} " + keyWord + " pairs, ";
	      switch(ring){
	      case 0:
		title += "Top ring";
		break;
	      case 1:
		title += "Middle ring";
		break;
	      case 2:
		title += "Bottom ring";
		break;
	      }
	      title += TString::Format(", fit iteration %d", fitInd);
	      title += "; #phi_{wave} (Radians); #Delta#deltat (ns)";
	      gr->SetTitle(title);
	      gr->GetXaxis()->SetLimits(0, TMath::TwoPi());
	      gr->GetYaxis()->SetRangeUser(-1, 1);
	      p->SetTitle(title);
	      p->GetYaxis()->SetRangeUser(-0.5, 0.5);
	      // p->GetYaxis()->SetRangeUser(-0.2, 0.2);	    
	    
	    }
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
  TGraph* grAdjMean_2 = (TGraph*) f->Get("grAdjMean_2");
  TGraph* grVertMean_2 = (TGraph*) f->Get("grVertMean_2");
  TGraph* grAdjRms_2 = (TGraph*) f->Get("grAdjRms_2");
  TGraph* grVertRms_2 = (TGraph*) f->Get("grVertRms_2");
  TGraph* grAdjSlope_2 = (TGraph*) f->Get("grAdjSlope_2");
  TGraph* grVertSlope_2 = (TGraph*) f->Get("grVertSlope_2");
  TGraph* grAdjMean_3 = (TGraph*) f->Get("grAdjMean_3");
  TGraph* grVertMean_3 = (TGraph*) f->Get("grVertMean_3");
  TGraph* grAdjRms_3 = (TGraph*) f->Get("grAdjRms_3");
  TGraph* grVertRms_3 = (TGraph*) f->Get("grVertRms_3");
  TGraph* grAdjSlope_3 = (TGraph*) f->Get("grAdjSlope_3");
  TGraph* grVertSlope_3 = (TGraph*) f->Get("grVertSlope_3");
  TGraph* grAdjMean_4 = (TGraph*) f->Get("grAdjMean_4");
  TGraph* grVertMean_4 = (TGraph*) f->Get("grVertMean_4");
  TGraph* grAdjRms_4 = (TGraph*) f->Get("grAdjRms_4");
  TGraph* grVertRms_4 = (TGraph*) f->Get("grVertRms_4");
  TGraph* grAdjSlope_4 = (TGraph*) f->Get("grAdjSlope_4");
  TGraph* grVertSlope_4 = (TGraph*) f->Get("grVertSlope_4");


  
  auto c1 = new TCanvas("fitCan", "fitCan", 1600, 1600);
  c1->Divide(2, 3);

  // TString iter0 = "Prior to fit";
  // TString iter1 = "Fitting t";
  // TString iter2 = "Fitting r, #phi";
  // TString iter3 = "Fitting t again";
  // TString iter4 = "Fitting z";    

  // TString iter0 = "Prior to fit";
  // TString iter1 = "Fitting #phi";

  TString iter0 = "Prior to fit";
  TString iter1 = "Fitting r, #phi";
  
  c1->cd(1);
  grAdjMean_0->SetTitle(iter0);
  grAdjMean_1->SetTitle(iter1);
  // grAdjMean_2->SetTitle(iter2);
  // grAdjMean_3->SetTitle(iter3);
  // grAdjMean_4->SetTitle(iter4);

  grAdjMean_0->SetMarkerColor(kBlue);
  grAdjMean_1->SetMarkerColor(kRed);
  // grAdjMean_2->SetMarkerColor(kBlack);
  // grAdjMean_3->SetMarkerColor(kMagenta);
  // grAdjMean_4->SetMarkerColor(kCyan);  

  grAdjMean_0->SetMarkerStyle(2);
  grAdjMean_1->SetMarkerStyle(2);
  // grAdjMean_2->SetMarkerStyle(2);
  // grAdjMean_3->SetMarkerStyle(2);
  // grAdjMean_4->SetMarkerStyle(2);  

  grAdjMean_0->SetLineColor(0);
  grAdjMean_1->SetLineColor(0);
  // grAdjMean_2->SetLineColor(0);
  // grAdjMean_3->SetLineColor(0);
  // grAdjMean_4->SetLineColor(0);  
  
  grAdjMean_0->Draw("ap");
  grAdjMean_1->Draw("psame");
  // grAdjMean_2->Draw("psame");
  // grAdjMean_3->Draw("psame");
  // grAdjMean_4->Draw("psame");  

  auto l1_1 = gPad->BuildLegend(0.8, 0.6, 1, 1);

  grAdjMean_0->SetTitle("Adjacent pairs means at each stage of iteration; Pair index; Mean (ns)");
  grAdjMean_0->GetYaxis()->SetRangeUser(-0.4, 0.4);
  
  c1->cd(3);
  grAdjRms_0->SetTitle(iter0);
  grAdjRms_1->SetTitle(iter1);
  // grAdjRms_2->SetTitle(iter2);
  // grAdjRms_3->SetTitle(iter3);
  // grAdjRms_4->SetTitle(iter4);  

  grAdjRms_0->SetMarkerColor(kBlue);
  grAdjRms_1->SetMarkerColor(kRed);
  // grAdjRms_2->SetMarkerColor(kBlack);
  // grAdjRms_3->SetMarkerColor(kMagenta);
  // grAdjRms_4->SetMarkerColor(kCyan);  

  grAdjRms_0->SetMarkerStyle(2);
  grAdjRms_1->SetMarkerStyle(2);
  // grAdjRms_2->SetMarkerStyle(2);
  // grAdjRms_3->SetMarkerStyle(2);
  // grAdjRms_4->SetMarkerStyle(2);  

  grAdjRms_0->SetLineColor(0);
  grAdjRms_1->SetLineColor(0);
  // grAdjRms_2->SetLineColor(0);
  // grAdjRms_3->SetLineColor(0);
  // grAdjRms_4->SetLineColor(0);      

  grAdjRms_0->Draw("ap");   
  grAdjRms_1->Draw("psame");
  // grAdjRms_2->Draw("psame");
  // grAdjRms_3->Draw("psame");
  // grAdjRms_4->Draw("psame");   
  
  auto l1_3 = gPad->BuildLegend(0.8, 0.6, 1, 1);
  grAdjRms_0->SetTitle("Adjacent pairs RMS at each stage of iteration; Pair index; RMS (ns)");
  grAdjRms_0->GetYaxis()->SetRangeUser(0, 0.1);

  c1->cd(5);  
  grAdjSlope_0->SetTitle(iter0);
  grAdjSlope_1->SetTitle(iter1);
  // grAdjSlope_2->SetTitle(iter2);
  // grAdjSlope_3->SetTitle(iter3);
  // grAdjSlope_4->SetTitle(iter4);  

  grAdjSlope_0->SetMarkerColor(kBlue);
  grAdjSlope_1->SetMarkerColor(kRed);
  // grAdjSlope_2->SetMarkerColor(kBlack);
  // grAdjSlope_3->SetMarkerColor(kMagenta);
  // grAdjSlope_4->SetMarkerColor(kCyan);  

  grAdjSlope_0->SetMarkerStyle(2);
  grAdjSlope_1->SetMarkerStyle(2);
  // grAdjSlope_2->SetMarkerStyle(2);
  // grAdjSlope_3->SetMarkerStyle(2);
  // grAdjSlope_4->SetMarkerStyle(2);  

  grAdjSlope_0->SetLineColor(0);
  grAdjSlope_1->SetLineColor(0);
  // grAdjSlope_2->SetLineColor(0);
  // grAdjSlope_3->SetLineColor(0);
  // grAdjSlope_4->SetLineColor(0);      

  grAdjSlope_0->Draw("ap");   
  grAdjSlope_1->Draw("psame");
  // grAdjSlope_2->Draw("psame");
  // grAdjSlope_3->Draw("psame");
  // grAdjSlope_4->Draw("psame");   
  
  auto l1_5 = gPad->BuildLegend(0.8, 0.6, 1, 1);
  grAdjSlope_0->SetTitle("Adjacent pairs slope at each stage of iteration; Pair index; Gradient (ns/rad)");
  grAdjSlope_0->GetYaxis()->SetRangeUser(-0.2, 0.2);



  c1->cd(2);
  grVertMean_0->SetTitle(iter0);
  grVertMean_1->SetTitle(iter1);
  // grVertMean_2->SetTitle(iter2);
  // grVertMean_3->SetTitle(iter3);
  // grVertMean_4->SetTitle(iter4);  

  grVertMean_0->SetMarkerColor(kBlue);
  grVertMean_1->SetMarkerColor(kRed);
  // grVertMean_2->SetMarkerColor(kBlack);
  // grVertMean_3->SetMarkerColor(kMagenta);
  // grVertMean_4->SetMarkerColor(kCyan);  

  grVertMean_0->SetMarkerStyle(2);
  grVertMean_1->SetMarkerStyle(2);
  // grVertMean_2->SetMarkerStyle(2);
  // grVertMean_3->SetMarkerStyle(2);
  // grVertMean_4->SetMarkerStyle(2);  

  grVertMean_0->SetLineColor(0);
  grVertMean_1->SetLineColor(0);
  // grVertMean_2->SetLineColor(0);
  // grVertMean_3->SetLineColor(0);
  // grVertMean_4->SetLineColor(0);  
  
  grVertMean_0->Draw("ap");
  grVertMean_1->Draw("psame");
  // grVertMean_2->Draw("psame");
  // grVertMean_3->Draw("psame");
  // grVertMean_4->Draw("psame");  

  auto l1_2 = gPad->BuildLegend(0.8, 0.6, 1, 1);
  grVertMean_0->SetTitle("Vertical pairs mean at each stage of iteration; Pair index; Mean (ns)");
  grVertMean_0->GetYaxis()->SetRangeUser(-0.4, 0.4);
  
  c1->cd(4);
  grVertRms_0->SetTitle(iter0);
  grVertRms_1->SetTitle(iter1);
  // grVertRms_2->SetTitle(iter2);
  // grVertRms_3->SetTitle(iter3);
  // grVertRms_4->SetTitle(iter4);  

  grVertRms_0->SetMarkerColor(kBlue);
  grVertRms_1->SetMarkerColor(kRed);
  // grVertRms_2->SetMarkerColor(kBlack);
  // grVertRms_3->SetMarkerColor(kMagenta);
  // grVertRms_4->SetMarkerColor(kCyan);  

  grVertRms_0->SetMarkerStyle(2);
  grVertRms_1->SetMarkerStyle(2);
  // grVertRms_2->SetMarkerStyle(2);
  // grVertRms_3->SetMarkerStyle(2);
  // grVertRms_4->SetMarkerStyle(2);  

  grVertRms_0->SetLineColor(0);
  grVertRms_1->SetLineColor(0);
  // grVertRms_2->SetLineColor(0);
  // grVertRms_3->SetLineColor(0);
  // grVertRms_4->SetLineColor(0);      

  grVertRms_0->Draw("ap");   
  grVertRms_1->Draw("psame");
  // grVertRms_2->Draw("psame");
  // grVertRms_3->Draw("psame");
  // grVertRms_4->Draw("psame");   
  
  auto l1_4 = gPad->BuildLegend(0.8, 0.6, 1, 1);
  grVertRms_0->SetTitle("Vertical pairs RMS at each stage of iteration; Pair index; RMS (ns)");
  grVertRms_0->GetYaxis()->SetRangeUser(0, 0.1);


  c1->cd(6);  
  grVertSlope_0->SetTitle(iter0);
  grVertSlope_1->SetTitle(iter1);
  // grVertSlope_2->SetTitle(iter2);
  // grVertSlope_3->SetTitle(iter3);
  // grVertSlope_4->SetTitle(iter4);  

  grVertSlope_0->SetMarkerColor(kBlue);
  grVertSlope_1->SetMarkerColor(kRed);
  // grVertSlope_2->SetMarkerColor(kBlack);
  // grVertSlope_3->SetMarkerColor(kMagenta);
  // grVertSlope_4->SetMarkerColor(kCyan);  

  grVertSlope_0->SetMarkerStyle(2);
  grVertSlope_1->SetMarkerStyle(2);
  // grVertSlope_2->SetMarkerStyle(2);
  // grVertSlope_3->SetMarkerStyle(2);
  // grVertSlope_4->SetMarkerStyle(2);

  grVertSlope_0->SetLineColor(0);
  grVertSlope_1->SetLineColor(0);
  // grVertSlope_2->SetLineColor(0);
  // grVertSlope_3->SetLineColor(0);
  // grVertSlope_4->SetLineColor(0);

  grVertSlope_0->Draw("ap");   
  grVertSlope_1->Draw("psame");
  // grVertSlope_2->Draw("psame");
  // grVertSlope_3->Draw("psame");
  // grVertSlope_4->Draw("psame");   
  
  auto l1_6 = gPad->BuildLegend(0.8, 0.6, 1, 1);
  grVertSlope_0->SetTitle("Vertical pairs slope at each stage of iteration; Pair index; Gradient (ns/rad)");
  grVertSlope_0->GetYaxis()->SetRangeUser(-0.2, 0.2);  
  
  return;
  

}
