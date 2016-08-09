void drawLookForPitchRollOffsetPlots2(){



  auto pitchRollChain = new TChain("pitchRollTree");
  // pitchRollChain->Add("lookForPitchRollOffsetsPlots/lookForPitchRollOffsetsPlots_*");
  // pitchRollChain->Add("lookForPitchRollOffsetsPlotsCoarse_2016_01_28/lookForPitchRollOffsetsPlots_*");
  pitchRollChain->Add("lookForPitchRollOffsetsPlots_352_2016-01-29_14-42-30.root");  

  const Int_t numPitches = 11;
  const Double_t startPitch = -2;
  const Double_t deltaPitch = 0.4;
  Double_t pitchOffsets[numPitches] = {0};
  for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
    pitchOffsets[pitchInd] = startPitch + deltaPitch*pitchInd;
  }
  
  const Int_t numRolls = 11;
  const Double_t startRoll = -2;
  const Double_t deltaRoll = 0.4;
  Double_t rollOffsets[numRolls] = {0};
  for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
    rollOffsets[rollInd] = startRoll + deltaRoll*rollInd;
  }  

  const Int_t numHeadings = 11;
  const Double_t startHeading = -2; //-0.5; //-2;
  const Double_t deltaHeading = 0.4; //0.1; //0.4;
  Double_t headingOffsets[numRolls] = {0};
  for(Int_t headingInd=0; headingInd < numHeadings; headingInd++){
    headingOffsets[headingInd] = startHeading + deltaHeading*headingInd;
  }  

  // const Int_t numPitches = 11;
  // const Double_t startPitch = -0.1; //-2;
  // const Double_t deltaPitch = 0.1; //0.4;
  // Double_t pitchOffsets[numPitches] = {0};
  // for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
  //   pitchOffsets[pitchInd] = startPitch + deltaPitch*pitchInd;
  // }
  
  // const Int_t numRolls = 11;
  // const Double_t startRoll = -0.5; //-2;
  // const Double_t deltaRoll = 0.1; //0.4;
  // Double_t rollOffsets[numRolls] = {0};
  // for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
  //   rollOffsets[rollInd] = startRoll + deltaRoll*rollInd;
  // }  


  std::vector<TProfile*> profs;  
  for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
    for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
      const Int_t headingInd = 5;
      
      if(!((pitchInd==5 && rollInd==5) || (pitchInd==6 && rollInd==5))){
      	continue;
      }
      
      TString profName = TString::Format("hProf_%d_%d", pitchInd, rollInd);
      TString profTitle = TString::Format("Pitch %2.1lf, Roll %2.1lf, Heading %2.1lf",
  					  pitchOffsets[pitchInd], rollOffsets[rollInd],
					  headingOffsets[headingInd]);

      TProfile* hProf = new TProfile(profName, profTitle, 256, 0, 360);

      // TString command = TString::Format("deltaThetaDegs[%d][%d]:zoomPhiDeg>>", pitchInd, rollInd);
      // command += profName;
      TString command = TString::Format("deltaThetaDegs[%d][%d][%d]:zoomPhiDeg>>", headingInd, pitchInd, rollInd);
      command += profName;

      std::cout << headingInd << "\t" << pitchInd << "\t" << rollInd << "\t" << command.Data() << std::endl;

      TString cuts = TString::Format("TMath::Abs(deltaThetaDegs[%d][%d][%d]) < 10", headingInd, pitchInd, rollInd);
      pitchRollChain->Draw(command, cuts, "goff");

      profs.push_back(hProf);

      hProf->SetFillColor(0);
      
      // TString drawOpt = "same";
      // if(pitchInd==0 && rollInd==0){
      // 	drawOpt = "";
      // }
      // hProf->Draw(drawOpt);
    }
  }

  TCanvas* c1 = RootTools::drawArrayOfHistosPrettily((TH1D**) &profs[0], (Int_t) profs.size(), NULL, NULL);
  TLegend* l = c1->BuildLegend();
    
  l->Draw();
  profs[0]->SetTitle("Pitch roll offsets");
  return 0;



  // for(int headingInd=0; headingInd < numHeadings; headingInd++){
  const Int_t startHeadingInd = 4;
  const Int_t endHeadingInd = 6;  

  ProgressBar p(numPitches*numRolls*(endHeadingInd-startHeadingInd));

  for(int headingInd=startHeadingInd; headingInd <= endHeadingInd; headingInd++){    
    TString title = TString::Format("RMS of #delta#phi distributions (Heading = %2.1lf); Pitch offset (degrees); Roll offset (degrees); RMS (Degrees)",headingOffsets[headingInd]);

    TH2D* hPhiRmsVals = new TH2D(TString::Format("hPhiRmsVals_%d", headingInd),
				 title,
				 numPitches, pitchOffsets[0], pitchOffsets[numPitches-1]+deltaPitch,
				 numRolls, rollOffsets[0], rollOffsets[numRolls-1]+deltaRoll);

    for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
      for(Int_t rollInd=0; rollInd < numRolls; rollInd++){

	// if(!((pitchInd==5 && rollInd==5) || (pitchInd==10 && rollInd==10))){
	// 	continue;
	// }

      
	TString histName = TString::Format("h_%d_%d_%d", headingInd, pitchInd, rollInd);	
	TString histTitle = TString::Format("Pitch %2.1lf, Roll %2.1lf, Heading %2.1lf",
					    pitchOffsets[pitchInd], rollOffsets[rollInd],
					    headingOffsets[headingInd]);

	TH1D* h = new TH1D(histName, histTitle, 200, -10, 10);

	// TString command = TString::Format("deltaThetaDegs[%d][%d]:zoomPhiDeg>>", pitchInd, rollInd);
	// command += histName;
	TString command = TString::Format("deltaThetaDegs[%d][%d][%d]>>", headingInd, pitchInd, rollInd);
	command += histName;

	TString cuts = TString::Format("TMath::Abs(deltaThetaDegs[%d][%d][%d]) < 10", headingInd, pitchInd, rollInd);
	pitchRollChain->Draw(command, cuts, "goff");

	Double_t rms = h->GetRMS();
	Double_t mean = h->GetMean();
	Double_t var = rms*rms;
	Double_t chiSq = var + mean*mean;
	hPhiRmsVals->SetBinContent(pitchInd+1, rollInd+1, chiSq);
	p++;
      }
    }
    new TCanvas();
    hPhiRmsVals->Draw("colz");
  }
  

}
