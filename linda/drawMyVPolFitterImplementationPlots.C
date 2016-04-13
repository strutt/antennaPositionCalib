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
  for(int i=0; i<n; i++){
    Double_t x = -xMax + dx*i;
    Double_t y = (10*x*x + rndm.Gaus(0, 1))/norm;
    gr->SetPoint(i, x, y);
  }
  
  return gr;
}

void drawMyVPolFitterImplementationPlots(){

  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-17-57.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_11-51-19.root");
  // TFile* f = TFile::Open("myVPolFitterImplementationPlots_2016-04-13_12-08-13.root");

  // TString fName = "myVPolFitterImplementationPlots_2016-04-13_14-28-17.root";
  // TString fName = "myVPolFitterImplementationPlots_2016-04-13_15-10-19.root";
  TString fName = "myVPolFitterImplementationPlots_2016-04-13_17-52-53.root";

  TFile* f = TFile::Open(fName);

  EColor cols[NUM_PHI] = {kRed, kBlue, kMagenta, kCyan, kGreen, kOrange, kViolet, kBlack,
			  kRed, kBlue, kMagenta, kCyan, kGreen, kOrange, kViolet, kBlack}; 

  const int numFitSteps = 5;
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
	    
	  }
	}
      }
    }
  }
  return;
  
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
    

  // auto c0 = new TCanvas();
  // grTotal_0->Draw();
  // grTotal->SetTitle("Sum of weighted figures of merit: Iteration 0");
  // grTotal_0->SetLineColor(kBlack);  
  // grAdjMean_0_weighted->Draw("lsame");
  // grAdjMean_0_weighted->SetTitle("Adjacent Mean: Iteration 0");
  // grAdjMean_0_weighted->SetLineColor(kBlue);
  // grVertMean_0_weighted->Draw("lsame");
  // grVertMean_0_weighted->SetLineColor(kRed);
  // grVertMean_0_weighted->SetTitle("Vertical Mean: Iteration 0");  
  // grAdjRms_0_weighted->Draw("lsame");
  // grAdjRms_0_weighted->SetLineColor(kMagenta);
  // grAdjRms_0_weighted->SetTitle("Adjacent RMS: Iteration 0");
  // grVertRms_0_weighted->Draw("lsame");
  // grVertRms_0_weighted->SetLineColor(kOrange);
  // grAdjSlope_0_weighted->Draw("lsame");
  // grAdjSlope_0_weighted->SetLineColor(kCyan);
  // grVertSlope_0_weighted->Draw("lsame");
  // grVertSlope_0_weighted->SetLineColor(kSpring);
  // c0->BuildLegend();
  // c0->SetLogy(1);
  
  // auto c01 = new TCanvas();
  // grTotal_1->Draw();
  // grTotal_1->SetLineColor(kBlack);  
  // grAdjMean_1_weighted->Draw("lsame");
  // grAdjMean_1_weighted->SetLineColor(kBlue);
  // grVertMean_1_weighted->Draw("lsame");
  // grVertMean_1_weighted->SetLineColor(kRed);
  // grAdjRms_1_weighted->Draw("lsame");
  // grAdjRms_1_weighted->SetLineColor(kMagenta);
  // grVertRms_1_weighted->Draw("lsame");
  // grVertRms_1_weighted->SetLineColor(kOrange);
  // grAdjSlope_1_weighted->Draw("lsame");
  // grAdjSlope_1_weighted->SetLineColor(kCyan);
  // grVertSlope_1_weighted->Draw("lsame");
  // grVertSlope_1_weighted->SetLineColor(kSpring);
  // c01->BuildLegend();
  // c01->SetLogy(1);


  auto c1 = new TCanvas();

  grAdjMean_0->Draw("ap");
  grAdjMean_1->Draw("psame");
  grAdjMean_2->Draw("psame");
  grAdjMean_3->Draw("psame");
  grVertMean_0->Draw("psame");  
  grVertMean_1->Draw("psame");  
  grVertMean_2->Draw("psame");  
  grVertMean_3->Draw("psame");  

  grAdjMean_0->SetMarkerColor(cols[0]);
  grAdjMean_1->SetMarkerColor(cols[1]);
  grAdjMean_2->SetMarkerColor(cols[2]);
  grAdjMean_3->SetMarkerColor(cols[3]);  
  grVertMean_0->SetMarkerColor(cols[4]);
  grVertMean_1->SetMarkerColor(cols[5]);
  grVertMean_2->SetMarkerColor(cols[6]);
  grVertMean_3->SetMarkerColor(cols[7]);  

  grAdjMean_0->SetLineColor(0);
  grAdjMean_1->SetLineColor(0);
  grAdjMean_2->SetLineColor(0);
  grAdjMean_3->SetLineColor(0);
  grVertMean_0->SetLineColor(0);
  grVertMean_1->SetLineColor(0);
  grVertMean_2->SetLineColor(0);
  grVertMean_3->SetLineColor(0);
  
  grAdjMean_0->SetMarkerStyle(2);
  grAdjMean_1->SetMarkerStyle(2);
  grVertMean_0->SetMarkerStyle(2);
  grVertMean_1->SetMarkerStyle(2);
  grAdjMean_2->SetMarkerStyle(2);
  grAdjMean_3->SetMarkerStyle(2);
  grVertMean_2->SetMarkerStyle(2);
  grVertMean_3->SetMarkerStyle(2);  

  grAdjMean_0->SetTitle("Adjacent Pairs: Iteration 0");
  grVertMean_0->SetTitle("Vertical Pairs: Iteration 0");
  grAdjMean_1->SetTitle("Adjacent Pairs: Iteration 1");  
  grVertMean_1->SetTitle("Vertical Pairs: Iteration 1");  
  grAdjMean_2->SetTitle("Adjacent Pairs: Iteration 2");
  grVertMean_2->SetTitle("Vertical Pairs: Iteration 2");
  grAdjMean_3->SetTitle("Adjacent Pairs: Iteration 3");  
  grVertMean_3->SetTitle("Vertical Pairs: Iteration 3");  


  c1->BuildLegend();
  grAdjMean_0->SetTitle("Mean #deltat_{measured} - #deltat_{expected} (ns); Antenna; #Delta#deltat (ns)");
  grAdjMean_0->GetYaxis()->SetRangeUser(-0.5, 0.5);
  c1->Update();

  auto c2 = new TCanvas();
  grAdjRms_0->Draw("ap");
  grAdjRms_1->Draw("psame");
  grAdjRms_2->Draw("psame");
  grAdjRms_3->Draw("psame");
  grVertRms_0->Draw("psame");  
  grVertRms_1->Draw("psame");  
  grVertRms_2->Draw("psame");  
  grVertRms_3->Draw("psame");  

  grAdjRms_0->SetMarkerColor(cols[0]);
  grAdjRms_1->SetMarkerColor(cols[1]);
  grAdjRms_2->SetMarkerColor(cols[2]);
  grAdjRms_3->SetMarkerColor(cols[3]);  
  grVertRms_0->SetMarkerColor(cols[4]);
  grVertRms_1->SetMarkerColor(cols[5]);
  grVertRms_2->SetMarkerColor(cols[6]);
  grVertRms_3->SetMarkerColor(cols[7]);  

  grAdjRms_0->SetLineColor(0);
  grAdjRms_1->SetLineColor(0);
  grAdjRms_2->SetLineColor(0);
  grAdjRms_3->SetLineColor(0);
  grVertRms_0->SetLineColor(0);
  grVertRms_1->SetLineColor(0);
  grVertRms_2->SetLineColor(0);
  grVertRms_3->SetLineColor(0);
  
  grAdjRms_0->SetMarkerStyle(2);
  grAdjRms_1->SetMarkerStyle(2);
  grVertRms_0->SetMarkerStyle(2);
  grVertRms_1->SetMarkerStyle(2);
  grAdjRms_2->SetMarkerStyle(2);
  grAdjRms_3->SetMarkerStyle(2);
  grVertRms_2->SetMarkerStyle(2);
  grVertRms_3->SetMarkerStyle(2);  
  grAdjRms_0->SetTitle("Adjacent Pairs: Iteration 0");
  grVertRms_0->SetTitle("Vertical Pairs: Iteration 0");
  grAdjRms_1->SetTitle("Adjacent Pairs: Iteration 1");  
  grVertRms_1->SetTitle("Vertical Pairs: Iteration 1");  
  grAdjRms_2->SetTitle("Adjacent Pairs: Iteration 2");
  grVertRms_2->SetTitle("Vertical Pairs: Iteration 2");
  grAdjRms_3->SetTitle("Adjacent Pairs: Iteration 3");  
  grVertRms_3->SetTitle("Vertical Pairs: Iteration 3");  


  c2->BuildLegend();
  grAdjRms_0->SetTitle("RMS #deltat_{measured} - #deltat_{expected} (ns); Antenna; RMS(#Delta#deltat) (ns)");
  grAdjRms_0->GetYaxis()->SetRangeUser(0, 0.5);
  c2->Update();   

  auto c3 = new TCanvas();

  grAdjSlope_0->Draw("ap");
  grAdjSlope_1->Draw("psame");
  grAdjSlope_2->Draw("psame");
  grAdjSlope_3->Draw("psame");
  grVertSlope_0->Draw("psame");  
  grVertSlope_1->Draw("psame");  
  grVertSlope_2->Draw("psame");  
  grVertSlope_3->Draw("psame");  

  grAdjSlope_0->SetMarkerColor(cols[0]);
  grAdjSlope_1->SetMarkerColor(cols[1]);
  grAdjSlope_2->SetMarkerColor(cols[2]);
  grAdjSlope_3->SetMarkerColor(cols[3]);  
  grVertSlope_0->SetMarkerColor(cols[4]);
  grVertSlope_1->SetMarkerColor(cols[5]);
  grVertSlope_2->SetMarkerColor(cols[6]);
  grVertSlope_3->SetMarkerColor(cols[7]);  

  grAdjSlope_0->SetLineColor(0);
  grAdjSlope_1->SetLineColor(0);
  grAdjSlope_2->SetLineColor(0);
  grAdjSlope_3->SetLineColor(0);
  grVertSlope_0->SetLineColor(0);
  grVertSlope_1->SetLineColor(0);
  grVertSlope_2->SetLineColor(0);
  grVertSlope_3->SetLineColor(0);
  
  grAdjSlope_0->SetMarkerStyle(2);
  grAdjSlope_1->SetMarkerStyle(2);
  grVertSlope_0->SetMarkerStyle(2);
  grVertSlope_1->SetMarkerStyle(2);
  grAdjSlope_2->SetMarkerStyle(2);
  grAdjSlope_3->SetMarkerStyle(2);
  grVertSlope_2->SetMarkerStyle(2);
  grVertSlope_3->SetMarkerStyle(2);  
  
  grAdjSlope_0->SetTitle("Adjacent Pairs: Iteration 0");
  grVertSlope_0->SetTitle("Vertical Pairs: Iteration 0");
  grAdjSlope_1->SetTitle("Adjacent Pairs: Iteration 1");  
  grVertSlope_1->SetTitle("Vertical Pairs: Iteration 1");  
  grAdjSlope_2->SetTitle("Adjacent Pairs: Iteration 2");
  grVertSlope_2->SetTitle("Vertical Pairs: Iteration 2");
  grAdjSlope_3->SetTitle("Adjacent Pairs: Iteration 3");  
  grVertSlope_3->SetTitle("Vertical Pairs: Iteration 3");  

  c3->BuildLegend();
  grAdjSlope_0->SetTitle("Slope #deltat_{measured} - #deltat_{expected} (ns); Antenna; Gradient(#Delta#deltat) (ns/rad)");
  grAdjSlope_0->GetYaxis()->SetRangeUser(-0.5, 0.5);
  c3->Update();
}    

  

