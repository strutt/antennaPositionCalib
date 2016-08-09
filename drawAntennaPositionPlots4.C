void drawAntennaPositionPlots4(){


  const int numFiles = 2;
  // TString lindaFiles[numFiles] = {"newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV2_2016_01_21_time_14_53_41.txt", "newLindaNumbers_LDBHPOL_NEW10_cosminV2_2016_01_21_time_14_03_25.txt"};				  
  TString lindaFiles[numFiles] = {"newLindaNumbers_4steps_WAISHPOL_NEW10_2016_01_19_time_15_02_11.txt",
				  "newLindaNumbers_4steps_WAISHPOL_NEW10_cosminV3_2016_01_25_time_12_24_25.txt"};
  TString shortLegNames[numFiles] = {"19th Jan", "cosminV3 25th Jan"};
  
  std::vector<Double_t> phiDegs[numFiles];
  auto c1 = new TCanvas();
  TGraph* grTop = NULL;
  for(int fileInd=0; fileInd < numFiles; fileInd++){
    // new TCanvas();
    auto pol = AnitaPol::kVertical;
    CrossCorrelator::directlyInsertGeometry(lindaFiles[fileInd], pol);
    auto geom = AnitaGeomTool::Instance();

    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      phiDegs[fileInd].push_back(TMath::RadToDeg()*geom->getAntPhiPositionRelToAftFore(ant, pol));
    }
    
    TGraph* grRings[AnitaRing::kNotARing];
    for(int ring=AnitaRing::kNotARing-1; ring>=0; ring--){
      grRings[ring] = new TGraph();
      switch(ring){
      case 0:
	grRings[ring]->SetTitle(shortLegNames[fileInd]+" Top Ring");
	break;
      case 1:
	grRings[ring]->SetTitle(shortLegNames[fileInd]+" Middle Ring");
	break;
      case 2:
	grRings[ring]->SetTitle(shortLegNames[fileInd]+" Bottom Ring");
	break;
      }
	
      for(int phi=0; phi<NUM_PHI; phi++){
	int ant = phi + ring*NUM_PHI;

	Double_t phiRad = geom->getAntPhiPositionRelToAftFore(ant, pol);
	Double_t r = geom->getAntR(ant, pol);    

	Double_t y = r*sin(phiRad);
	Double_t x = r*cos(phiRad);    

	grRings[ring]->SetPoint(grRings[ring]->GetN(), x, y);
      }

      grRings[ring]->SetMarkerStyle(fileInd+2);
      grRings[ring]->SetMarkerColor(1+ring);
      grRings[ring]->SetFillColor(0);
      grRings[ring]->SetLineColor(0);                
      TString opt = ring==AnitaRing::kNotARing-1 && fileInd==0 ? "ap" : "psame";
      if(ring==AnitaRing::kNotARing-1 && fileInd==0){
	grTop = grRings[ring];
      }
      grRings[ring]->Draw(opt);
    }
  }
  auto l = c1->BuildLegend();
  l->Draw();
  grTop->SetTitle("Antenna phase center positions from "+shortLegNames[0]+" and "+shortLegNames[1]+" fits; x (m); y (m)");
  
  TGraph* grDeltaPhiDegs = new TGraph();
  Double_t sum=0;
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    double dphi = RootTools::getDeltaAngleDeg(phiDegs[1].at(ant), phiDegs[0].at(ant));
    grDeltaPhiDegs->SetPoint(ant, ant, dphi);
    sum += dphi;
  }

  cout << sum / NUM_SEAVEYS << endl;
  new TCanvas();
  grDeltaPhiDegs->Draw("alp");
  grDeltaPhiDegs->SetTitle("#phi position "+shortLegNames[1]+ " - "+shortLegNames[0]+" phase centers; Antenna; #delta#phi (Degrees)");  
}
