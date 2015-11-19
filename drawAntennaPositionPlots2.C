{
  // const Int_t NUM_SEAVEYS = 48;
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

  std::vector<Double_t> ants;
  std::vector<Double_t> photoGrammetryRs;
  std::vector<Double_t> myFittedRs;
  std::vector<Double_t> deltaRs;

  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    ants.push_back(ant);

    geom->useKurtAnitaIIINumbers(1);
    photoGrammetryRs.push_back(geom->getAntR(ant));

    geom->useKurtAnitaIIINumbers(0);
    myFittedRs.push_back(geom->getAntR(ant));

    deltaRs.push_back(fittedDeltaRs[ant]);
  }

  auto gr0 = new TGraph(NUM_SEAVEYS, &ants[0], &myFittedRs[0]);
  auto gr1 = new TGraph(NUM_SEAVEYS, &ants[0], &photoGrammetryRs[0]);

  gr0->SetTitle("Antenna radial shift; Antenna Number; Radius (m)");

  TGraph* grs[2] = {gr0, gr1};
  TString titles0[2] = {"Photogrammetry", "Fitted phase center"};
  TLegend* l0 = RootTools::makeLegend(grs, 2, titles0, "p");
  
  
  auto c = new TCanvas();
  gr0->SetMarkerColor(kRed);
  gr0->SetMarkerStyle(8);
  gr0->Draw("ap");
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(8);
  gr1->Draw("psame");
  l0->Draw();

  
  auto c2 = new TCanvas();
  auto gr2 = new TGraph(NUM_SEAVEYS, &ants[0], &deltaRs[0]);
  gr2->SetMarkerStyle(8);
  gr2->Draw("ap");
  gr2->GetYaxis()->SetNoExponent(1);
  gr2->SetTitle("Antenna radial shift; Antenna Number; #deltaR (m)");

  const int NUM_RINGS = 3;
  TString titles[NUM_RINGS] = {"Top ring", "Middle ring", "Bottom ring"};
  auto c3 = new TCanvas();
  TGraph* gr2s[NUM_RINGS];
  for(int ring=0; ring<NUM_RINGS; ring++){
    gr2s[ring] = new TGraph(NUM_PHI, &ants[0], &deltaRs[ring*NUM_PHI]);
    gr2s[ring]->SetMarkerStyle(8);
    gr2s[ring]->SetTitle("Antenna radial shift; Phi-sector; #deltaR (m)");
  }
  RootTools::drawArrayOfTGraphsPrettily(gr2s, NUM_RINGS, "p", c3);
  TLegend* l = RootTools::makeLegend(gr2s, NUM_RINGS, titles, "p");
  l->Draw();
  gr2s[0]->SetMaximum(-0.09);
  gr2s[0]->SetMinimum(-0.1);
  gr2->SetMaximum(-0.09);
  gr2->SetMinimum(-0.1);
  
}
