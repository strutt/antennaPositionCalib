void plotConvergence(int step = 1){

  int count = 0;
  int tempcount=0;

  ifstream myInput(Form("converge_steps_%i.txt", step));

//    ifstream myInput("lindaMacros/ANITAsymmetric.txt");

  double meanH, rmsH, gradH, meanV, rmsV, gradV;

  double sumMeanH[1000];
  double sumRmsH[1000];
  double sumGradH[1000];
  double sumMeanV[1000];
  double sumRmsV[1000];
  double sumGradV[1000];
  double x[1000];

  if (myInput.is_open()){
    while(myInput.good()){
      
      myInput >> meanH >> rmsH >> gradH >> meanV >> rmsV >> gradV;
      
      if (tempcount%100==0){
	sumMeanH[count] = meanH;
	sumRmsH[count] = rmsH;
	sumGradH[count] = gradH;
	sumMeanV[count] = meanV;
	sumRmsV[count] = rmsV;
	sumGradV[count] = gradV;
	x[count] = tempcount*1.;
	count++;
	cout << meanH << " " << rmsH << " " << gradH << endl;
      }
      
      tempcount++;
    }
  }
  
  
  count--;

  TGraph *gMeanH = new TGraph(count, x, sumMeanH);
  TGraph *gRmsH  = new TGraph(count, x, sumRmsH);
  TGraph *gGradH = new TGraph(count, x, sumGradH);
  TGraph *gMeanV = new TGraph(count, x, sumMeanV);
  TGraph *gRmsV  = new TGraph(count, x, sumRmsV);
  TGraph *gGradV = new TGraph(count, x, sumGradV);

  gMeanH->SetLineColor(kCyan+1);
  gMeanH->SetLineWidth(2);

  gRmsH->SetLineColor(kMagenta+3);
  gRmsH->SetLineStyle(7);
  gRmsH->SetLineWidth(2);

  gGradH->SetLineColor(kGreen+2);
  gGradH->SetLineWidth(2);
  gGradH->SetLineStyle(8);


  gMeanV->SetLineColor(kCyan+1);
  gMeanV->SetLineWidth(2);

  gRmsV->SetLineColor(kMagenta+3);
  gRmsV->SetLineStyle(7);
  gRmsV->SetLineWidth(2);

  gGradV->SetLineColor(kGreen+2);
  gGradV->SetLineWidth(2);
  gGradV->SetLineStyle(8);


  TCanvas *c1 = new TCanvas("c1");
 

  gMeanH->Draw("Al");
  gMeanH->SetTitle("Convergence:horizontal antenna pairs;Iteration;Value");

  gMeanH->GetYaxis()->SetRangeUser(0.01,10000);

  c1->SetLogy();
  gMeanH->Draw("al");

  gRmsH->Draw("l");
  gGradH->Draw("l");

  TLegend *leg = new TLegend(0.89, 0.89, 0.55, 0.7);
  leg->SetFillColor(kWhite);
  leg->AddEntry(gMeanH, "#Sigma MEAN", "l");
  leg->AddEntry(gRmsH,  "#Sigma RMS" , "l");
  leg->AddEntry(gGradH, "#Sigma GRAD", "l");

  leg->Draw();

  c1->Print(Form("Convergence_HorizontalAntennaPairs_steps_%i.png", step));
  c1->Print(Form("Convergence_HorizontalAntennaPairs_steps_%i.pdf", step));


  gMeanV->Draw("Al");
  gMeanV->SetTitle("Convergence: vertical antenna pairs;Iteration;Value");

  gMeanV->GetYaxis()->SetRangeUser(0.01,10000);

  c1->SetLogy();
  gMeanV->Draw("al");

  gRmsV->Draw("l");
  gGradV->Draw("l");

  leg->Draw();

  c1->Print(Form("Convergence_VerticalAntennaPairs_steps_%i.png", step));
  c1->Print(Form("Convergence_VerticalAntennaPairs_steps_%i.pdf", step));


  
}
