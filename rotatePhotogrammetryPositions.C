void rotatePhotogrammetryPositions(){


  auto geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
  
  Double_t photoX[NUM_SEAVEYS];
  Double_t photoY[NUM_SEAVEYS];
  Double_t photoZ[NUM_SEAVEYS];

  Double_t rotX[NUM_SEAVEYS];
  Double_t rotY[NUM_SEAVEYS];
  Double_t rotZ[NUM_SEAVEYS];

  Double_t newAntR[NUM_SEAVEYS];
  Double_t newAntPhi[NUM_SEAVEYS];

  
  TGraph* grPhoto = new TGraph();
  TGraph* grRot = new TGraph();
  grPhoto->SetMarkerStyle(2);
  grRot->SetMarkerStyle(2);  
  grPhoto->SetMarkerColor(kBlue);
  grRot->SetMarkerColor(kRed);

  Double_t pitch = 1.6; //20; //1.6*10;
  Double_t roll = 2.0; //*10;
  Double_t heading = 0;

  ofstream bensFile("newBensNumbers_rotatingByHand.txt");
  
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
  // for(int ant=16; ant<32; ant++){    
    geom->getAntXYZ(ant, photoX[ant], photoY[ant], photoZ[ant], pol);

    // std::cout << photoX[ant] << "\t" << photoY[ant] << "\t" << photoZ[ant] << std::endl;



    TVector3 rollAxis=geom->fRollRotationAxis;
    TVector3 pitchAxis=geom->fPitchRotationAxis;

    TVector3 arbDir(photoX[ant], photoY[ant], photoZ[ant]);
     
    rollAxis.Rotate(heading*TMath::DegToRad(),geom->fHeadingRotationAxis);
    pitchAxis.Rotate(heading*TMath::DegToRad(),geom->fHeadingRotationAxis);
    rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     
    arbDir.Rotate(-1.*roll*TMath::DegToRad(),rollAxis);
    arbDir.Rotate(-1*pitch*TMath::DegToRad(),pitchAxis);

    // std::cout << roll << "\t" << pitch << "\t in function " << __PRETTY_FUNCTION__ << std::endl;
     
    arbDir.Rotate(-1*heading*TMath::DegToRad(),geom->fHeadingRotationAxis);

    //need to rotate roll and pitch axes - in this function heading has been
    //chosen as +ve x axis, pitch and roll axes are defined in terms of x and
    //y axes from ANITA setup document
    //rollAxis.Rotate(-135.*TMath::DegToRad(),geom->fHeadingRotationAxis);
    //pitchAxis.Rotate(45.*TMath::DegToRad(),geom->fHeadingRotationAxis);
    //pitchAxis.Rotate(-1*heading*TMath::DegToRad(),geom->fHeadingRotationAxis);

    rotX[ant] = arbDir.X();
    rotY[ant] = arbDir.Y();
    rotZ[ant] = arbDir.Z();


    grPhoto->SetPoint(grPhoto->GetN(), photoX[ant], photoY[ant]);
    grRot->SetPoint(grRot->GetN(), rotX[ant], rotY[ant]);


    // xAntFromVerticalHorn[ant]=rAntFromVerticalHorn[ant]*TMath::Cos(azCentreFromVerticalHorn[ant]);
    newAntR[ant] = TMath::Sqrt(rotX[ant]*rotX[ant] + rotY[ant]*rotY[ant]);
    newAntPhi[ant] = TMath::ATan2(rotY[ant], rotX[ant]);

    if(newAntPhi[ant] < 0){
      newAntPhi[ant] += TMath::TwoPi();
    }
    
    // std::cout << newAntR[ant] << "\t" << geom->getAntR(ant, pol) << std::endl;
    // // std::cout << newAntPhi[ant] << "\t" << geom->getAntPhiRelToAftFore(ant, pol) << std::endl;

    // // std::cout << "phis: " << newAntPhi[ant]*TMath::RadToDeg() << "\t" << geom->getAntPhiPositionRelToAftFore(ant, pol)*TMath::RadToDeg() << "\t"
    // // 	      << geom->getAntPhiPosition(ant, pol)*TMath::RadToDeg() << std::endl;    
    // std::cout << "phis: " << ant << "\t" << newAntPhi[ant]*TMath::RadToDeg() << "\t" 
    // 	      << geom->getAntPhiPosition(ant, pol)*TMath::RadToDeg() << std::endl;        


    Double_t deltaT = 0;
    Double_t deltaR = newAntR[ant] - geom->getAntR(ant, pol);
    Double_t deltaPhi = newAntPhi[ant] - geom->getAntPhiPosition(ant, pol);
    Double_t deltaZ = 0;    

    bensFile << ant << "  " << deltaR << "  " << deltaZ << "  " << deltaPhi << "  " << deltaT << std::endl;
  }

  
  grPhoto->Draw("ap");
  grRot->Draw("psame");  


  

}
