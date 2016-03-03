void drawLookForPitchRollOffsetPlots(){



  auto pitchRollChain = new TChain("pitchRollTree");
  pitchRollChain->Add("lookForPitchRollOffsetsPlots/lookForPitchRollOffsetsPlots_*");

  // const Int_t numPitches = 11;
  // const Double_t startPitch = -2;
  // const Double_t deltaPitch = 0.4;
  // Double_t pitchOffsets[numPitches] = {0};
  // for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
  //   pitchOffsets[pitchInd] = startPitch + deltaPitch*pitchInd;
  // }
  
  // const Int_t numRolls = 11;
  // const Double_t startRoll = -2;
  // const Double_t deltaRoll = 0.4;
  // Double_t rollOffsets[numRolls] = {0};
  // for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
  //   rollOffsets[rollInd] = startRoll + deltaRoll*rollInd;
  // }  


  const Int_t numPitches = 11;
  const Double_t startPitch = -0.1; //-2;
  const Double_t deltaPitch = 0.1; //0.4;
  Double_t pitchOffsets[numPitches] = {0};
  for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
    pitchOffsets[pitchInd] = startPitch + deltaPitch*pitchInd;
  }
  
  const Int_t numRolls = 11;
  const Double_t startRoll = -0.5; //-2;
  const Double_t deltaRoll = 0.1; //0.4;
  Double_t rollOffsets[numRolls] = {0};
  for(Int_t rollInd=0; rollInd < numRolls; rollInd++){
    rollOffsets[rollInd] = startRoll + deltaRoll*rollInd;
  }  


  // std::vector<TProfile*> profs;  
  // for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
  //   for(Int_t rollInd=0; rollInd < numRolls; rollInd++){

  //     if(!((pitchInd==5||pitchInd==6) && rollInd==5)){
  // 	continue;
  //     }
      
  //     TString profName = TString::Format("hProf_%d_%d", pitchInd, rollInd);
  //     TString profTitle = TString::Format("Pitch %2.1lf, Roll %2.1lf",
  // 					  pitchOffsets[pitchInd], rollOffsets[rollInd]);

  //     TProfile* hProf = new TProfile(profName, profTitle, 256, 0, 360);

  //     // TString command = TString::Format("deltaPhiDegs[%d][%d]:zoomPhiDeg>>", pitchInd, rollInd);
  //     // command += profName;
  //     TString command = TString::Format("deltaThetaDegs[%d][%d]:zoomPhiDeg>>", pitchInd, rollInd);
  //     command += profName;

  //     std::cout << pitchInd << "\t" << rollInd << "\t" << command.Data() << std::endl;

  //     TString cuts = TString::Format("TMath::Abs(deltaPhiDegs[%d][%d]) < 10", pitchInd, rollInd);
  //     pitchRollChain->Draw(command, cuts, "goff");

  //     profs.push_back(hProf);

  //     hProf->SetFillColor(0);
      
  //     // TString drawOpt = "same";
  //     // if(pitchInd==0 && rollInd==0){
  //     // 	drawOpt = "";
  //     // }
  //     // hProf->Draw(drawOpt);
  //   }
  // }

  // TCanvas* c1 = RootTools::drawArrayOfHistosPrettily((TH1D**) &profs[0], (Int_t) profs.size(), NULL, NULL);
  // TLegend* l = c1->BuildLegend();
    
  // l->Draw();
  // profs[0]->SetTitle("Pitch roll offsets");

  // return ;




  TH2D* hThetaRmsVals = new TH2D("hThetaRmsVals",
				 "RMS of #delta#theta distributions; Pitch offset (degrees); Roll offset (degrees); RMS (Degrees)",
				 numPitches, pitchOffsets[0], pitchOffsets[numPitches-1]+deltaPitch,
				 numRolls, rollOffsets[0], rollOffsets[numRolls-1]+deltaRoll);

  ProgressBar p(numPitches*numRolls);
  for(Int_t pitchInd=0; pitchInd < numPitches; pitchInd++){
    for(Int_t rollInd=0; rollInd < numRolls; rollInd++){

      TString histName = TString::Format("h_%d_%d", pitchInd, rollInd);
      TString histTitle = TString::Format("Pitch %2.1lf, Roll %2.1lf",
					  pitchOffsets[pitchInd], rollOffsets[rollInd]);

      TH1D* h = new TH1D(histName, histTitle, 200, -10, 10);

      // TString command = TString::Format("deltaPhiDegs[%d][%d]:zoomPhiDeg>>", pitchInd, rollInd);
      // command += histName;
      TString command = TString::Format("deltaThetaDegs[%d][%d]>>", pitchInd, rollInd);
      command += histName;

      TString cuts = TString::Format("TMath::Abs(deltaPhiDegs[%d][%d]) < 10", pitchInd, rollInd);
      pitchRollChain->Draw(command, cuts, "goff");


      hThetaRmsVals->SetBinContent(pitchInd+1, rollInd+1, h->GetRMS());
      p++;
    }
  }
  hThetaRmsVals->Draw("colz");
  

}
