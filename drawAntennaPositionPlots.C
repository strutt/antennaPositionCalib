{
  const Int_t NUM_SEAVEYS = 48;
  const Double_t extraCableDelays[NUM_SEAVEYS] =
    {0        , 0.139733  , 0.00697198, 0.0870043, 0.0502864, 0.139811 , 0.0564542, 0.143131 ,
     0.0755872, 0.10746   , -0.0529778, 0.109205 , 0.0199421, 0.062739 , 0.0813765, 0.033311 ,
     0.107176 , 0.141001  , 0.104494  , 0.16646  , 0.135572 , 0.14326  , 0.172411 , 0.113276 ,
     0.135898 , 0.0634696 , 0.078761  , 0.0890165, 0.198665 , 0.0649681, 0.057582 , 0.0800643,
     0.0635002, 0.211805  , 0.133592  , 0.136322 , 0.125028 , 0.0841643, 0.13517  , 0.0633426,
     0.0853349, 0.0491331 , 0.0996015 , 0.0681319, 0.0430019, 0.0380842, 0.0419707, -0.2944  };

  const Double_t fittedDeltaRs[NUM_SEAVEYS] =
    {-0.0959145, -0.0969593, -0.0976093, -0.0971427, -0.0972972, -0.0971141, -0.0971505, -0.0964771,
     -0.0960801, -0.094751 , -0.0938789, -0.0931123, -0.0941293, -0.0937031, -0.0944757, -0.0954,
     -0.0958407, -0.0967933, -0.0975012, -0.0970367, -0.0971045, -0.0971866, -0.0969005, -0.0965772,
     -0.0958794, -0.0949967, -0.0939146, -0.0936601, -0.0936972, -0.0940001, -0.0947252, -0.0954105,
     -0.0959116, -0.0966248, -0.0972106, -0.0971995, -0.0969414, -0.0970038, -0.0968214, -0.0965757,
     -0.0959988, -0.0949063, -0.0933371, -0.0934933, -0.0942268, -0.094005 , -0.0948142, -0.0957325};

  const Double_t fittedDeltaPhiDeg[NUM_SEAVEYS] =
    {0        , -1.60782  , -1.66883  , -1.67935  , -1.24366  , -0.693216  , -0.750744, -0.414814 ,
     -0.100636, -0.520359 , -0.45824  , -0.415373 , -0.610565 , -0.584207  , -0.912949, -1.27139  ,
     0.286926 , 0.440517  , 0.0755035 , -0.0280439, 0.0421618 , -0.0348197 , -0.225099, -0.187683 ,
     -0.106832, -0.0549485, -0.231417 , 0.0716926 , -0.217488 , -0.240067  , 0.240413 , 0.0476125 ,
     0.0949495, 0.185487  , -0.0496863, 0.112578  , -0.0353721, 9.89751e-06, -0.299654, -0.148599 ,
     -0.28825 , -0.32434  , -0.500945 , -0.0953967, -0.226163 , -0.312985  , 0.016517 , 0.00237409};

  const Double_t fittedDeltaZ[NUM_SEAVEYS] =
    {0          , -0.024893  , -0.0298469, -0.0213799, -0.0155314, -0.0220142, -0.0196334, -0.00451026,
     -0.00677221, -0.00123661, 0.011139  , 0.018916  , 0.0163255 , 0.0311236 , 0.0311453 , 0.0216939  ,
     0.0010381  , -0.029041  , -0.0388661, -0.0300287, -0.0284601, -0.0308645, -0.0232042, -0.00913271,
     -0.00738426, -0.00188341, 0.00968099, 0.0153015 , 0.0145462 , 0.0280903 , 0.0245277 , 0.0179969  ,
     0.00502674 , -0.0222594 , -0.0301759, -0.0255358, -0.0279802, -0.0314485, -0.0257633, -0.0153958 ,
     -0.0152704 , -0.0104062 , 0.00701353, 0.0128274 , 0.0100978 , 0.0285769 , 0.0277663 , 0.0182231  };

  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();

  const int numPhi = 16;
  const int numRings = 3;
  const int numRings2 = numRings + 1;

  TGraph* grs[numPhi];
  TGraph* grs0[numPhi];  

  Double_t means[numRings];
  for(int ring=0; ring<numRings2; ring++){
    means[ring] = 0;
  }
  
  for(int phi=0; phi<numPhi; phi++){
    grs[phi] = new TGraph();
    grs0[phi] = new TGraph();    
    for(int ring=0; ring<numRings; ring++){
      int ant = phi + numPhi*ring;
      grs[phi]->SetPoint(ring, geom->getAntPhiPosition(ant)*TMath::RadToDeg(), geom->getAntZ(ant));
      grs0[phi]->SetPoint(ring, geom->getAntPhiPosition(ant)*TMath::RadToDeg(), geom->getAntR(ant));      
      if(ring > 0){
	means[ring+1] += (geom->getAntZ(ant)/numPhi);
      }
      else{
	if((phi & 1)==0){
	  means[0] += (geom->getAntZ(ant)/(numPhi/2));
	}
	else{
	  means[1] += (geom->getAntZ(ant)/(numPhi/2));
	}
      }
    }
    // cout << grs[phi]->GetN() << endl;
    grs[phi]->SetMarkerStyle(8);
    grs[phi]->SetMarkerSize(2);
    grs0[phi]->SetMarkerStyle(8);
    grs0[phi]->SetMarkerSize(2);
  }


  TGraph* grs2[numRings2];
  TGraph* grs3[numRings2];  
  for(int ring=0; ring<numRings2; ring++){
    // means[ring]=0;
    grs2[ring] = new TGraph();
    grs2[ring]->SetMarkerStyle(8);
    grs2[ring]->SetMarkerSize(2);

    grs3[ring] = new TGraph();
    grs3[ring]->SetMarkerStyle(8);
    grs3[ring]->SetMarkerSize(2);
    
  }
  Int_t numPointsLower = 0;
  Int_t numPointsHigher = 0;  
  for(int ring=0; ring<numRings; ring++){
    for(int phi=0; phi<numPhi; phi++){
      if(ring > 0){
	grs2[ring+1]->SetPoint(phi, grs[phi]->GetX()[ring], grs[phi]->GetY()[ring]-means[ring+1]);

	Double_t dz = grs[phi]->GetY()[ring]-means[ring+1];
	Double_t r = grs0[phi]->GetY()[ring];
	Double_t theta = TMath::RadToDeg()*TMath::ATan2(dz, r);
	grs3[ring+1]->SetPoint(phi, grs[phi]->GetX()[ring], theta);

	cout << ring << "\t" << phi << "\t" << dz << "\t" << r << "\t" << theta << endl;
	
      }
      else{
	if((phi%2)==0){
	  Double_t dz = grs[phi]->GetY()[0]-means[0];
	  Double_t r = grs0[phi]->GetY()[0];
	  Double_t theta = TMath::RadToDeg()*TMath::ATan2(dz, r);

	  grs2[0]->SetPoint(numPointsHigher, grs[phi]->GetX()[0], dz);
	  grs3[0]->SetPoint(numPointsHigher, grs[phi]->GetX()[0], theta);
	  numPointsHigher++;
	}
	else{
	  Double_t dz = grs[phi]->GetY()[0]-means[1];
	  Double_t r = grs0[phi]->GetY()[0];
	  Double_t theta = TMath::RadToDeg()*TMath::ATan2(dz, r);
	
	  grs2[1]->SetPoint(numPointsLower, grs[phi]->GetX()[0], dz);
	  grs3[1]->SetPoint(numPointsLower, grs[phi]->GetX()[0], theta);
	  numPointsLower++;
	}
      }
    }
  }

  TCanvas* c1 = RootTools::drawArrayOfTGraphsPrettily(grs, numPhi, "p", NULL, NULL);
  grs[0]->GetXaxis()->SetLimits(-10, 370);
  grs[0]->SetTitle("Antenna height; Azimuth (degrees); Height (m)");
  grs[0]->GetYaxis()->SetNoExponent(1);  
  // for(int ring=0; ring<numRings2; ring++){
  //   for(int samp=0; samp < grs2[ring]->GetN(); samp++){
  //     cout << samp << "\t" << grs2[ring]->GetX()[samp] << "\t" << grs2[ring]->GetY()[samp] << endl;
  //   }
  // }
  

  TString titles1[numPhi];
  for(int phi=0; phi<numPhi; phi++){
    titles1[phi] = TString::Format("Phi-sector %d", phi);
  }
  TLegend* l1 = RootTools::makeLegend(grs, numPhi, titles1, "p", 0.85, 0.1, 1, 0.9);
  l1->Draw();
  

  
  TCanvas* c2 = RootTools::drawArrayOfTGraphsPrettily(grs2, numRings2, "lp", NULL, NULL);
  grs2[0]->GetXaxis()->SetLimits(-10, 370);
  grs2[0]->SetTitle("Antenna height relative to ring mean; Azimuth (degrees); Height (m)");  
  grs2[0]->GetYaxis()->SetNoExponent(1);
  TString titles2[numRings2] = {"Upper-top ring", "Lower-top ring", "Middle ring", "Bottom ring"};
  TLegend* l2 = RootTools::makeLegend(grs2, numRings2, titles2, "lp", 0.8, 0.8, 1, 1);
  l2->Draw();


  
  TCanvas* c3 = RootTools::drawArrayOfTGraphsPrettily(grs3, numRings2, "lp", NULL, NULL);
  grs3[0]->GetXaxis()->SetLimits(-10, 370);
  grs3[0]->SetTitle("Angle of antenna position relative to centre at mean ring height; Azimuth (degrees); Angle relative to mean z-axis position (Degrees)");
  grs3[0]->GetYaxis()->SetNoExponent(1);
  TString titles3[numRings2] = {"Upper-top ring", "Lower-top ring", "Middle ring", "Bottom ring"};
  TLegend* l3 = RootTools::makeLegend(grs3, numRings2, titles3, "lp", 0.8, 0.8, 1, 1);
  l3->Draw();
  


}
