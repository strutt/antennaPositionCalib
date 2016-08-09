{


  
  TFile *f = TFile::Open("makeCorrelationSummariesJustOneEventPlots_352_2016-02-25_18-08-09.root");
  for(int combo=0; combo<NUM_COMBOS; combo++){
    
    auto gr = (TGraph*) f->Get(TString::Format("gr_hZ%d", combo));
    if(gr){

      auto h = (TH2D*) f->Get(TString::Format("hZ%d", combo));
      if(gPad==NULL){
	h->Draw();
      }

      if(h){
	Double_t maxVal = h->GetMaximum();
	if(maxVal > 0.1){
	  gr->SetLineColor(gStyle->GetColorPalette(int(254*maxVal/0.5)));
	  if(gr->GetN() > 0){
	    gr->Draw("lsame");
	  }
	}
      }
    }
  }
}
