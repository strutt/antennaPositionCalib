{

  CrossCorrelator* cc = new CrossCorrelator();
  const int numFiles = 3;
  TString files[numFiles] = {"generateDeltaTLookupHistogramsPlotsLindaPhaseCenters.root",
			     "generateDeltaTLookupHistogramsPlotsPhotogrammetry.root",
			     "generateDeltaTLookupHistogramsPlots.root"};  
  TCanvas* c1[336] = {NULL};

  for(int fileInd=0; fileInd < numFiles; fileInd++){
    TFile* f = TFile::Open(files[fileInd]);

    Int_t counter=0;

    Double_t chi2 = 0;
    Double_t chi2Count = 0;
    for(int combo=0; combo<336; combo++){
      TGraph* gr = (TGraph*) f->Get(TString::Format("grCombo%d", combo));

      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      // if(ant2 - ant1 == 16 || ant2 - ant1 == 32){

      if( ((ant1/16) == (ant2/16) && (ant2 - ant1 == 1 || ant2 - ant1 == 15))
	  || (ant2 - ant1 == 16 || ant2 - ant1 == 32)
	 ){	 

	// std::cout << ant1 << "\t" << ant2 << std::endl;
	TString title = TString::Format("ant1 = %d, ant2 = %d", ant1, ant2);
      
	// if(combo >= 3) continue;
	TString fitName = TString::Format("f%d_%d_%d", ant1, ant2, fileInd);
	TF1* f1 = new TF1(fitName, "[0]+[1]*x");	
	if(fileInd==0){
	  gr->SetMarkerColor(kBlue);
	  f1->SetLineColor(kBlue);
	  title += ", Phase center numbers";
	}
	else if(fileInd==1){
	  gr->SetMarkerColor(kRed);
	  f1->SetLineColor(kRed);	  
	  title += ", Photogrammetry numbers";
	}
	else{
	  gr->SetMarkerColor(kBlack);
	  f1->SetLineColor(kBlack);
	  
	  title += ", Random";
	}
	gr->SetMarkerStyle(2);
	gr->SetTitle(title);


	for(int i=0; i<gr->GetN(); i++){
	  Double_t y = gr->GetY()[i];
	  chi2 += y*y;
	  chi2Count++;
	}

	if(counter==0){
	  if(c1[combo]==NULL){
	    c1[combo] = new TCanvas();
	    gr->Draw("ap");
	  }
	  else{
	    c1[combo]->cd();
	    gr->Draw("psame");
	  }

	  gr->GetXaxis()->SetLimits(0, 360);
	  // TF1* f1 = new TF1(gr->GetTitle(), "[0]");//+[1]*x");      
	  gr->Fit(f1, "Q");
	  // }
	  // else{
	  // 	gr->Draw("p");
	  // }
	}
	counter++;
      }
    }

    std::cout << f->GetTitle() << "\t" << TMath::Sqrt(chi2/chi2Count) << std::endl;
    f->Close();
  }

}
