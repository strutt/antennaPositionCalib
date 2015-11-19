{
   // show several TGaxis formats
  gROOT->Reset();

  // c1 = new TCanvas("c1","Examples of Gaxis",10,10,700,500);

  // c1->Range(-10,-1,10,1);

  // TGaxis *axis1 = new TGaxis(-5,-0.2,6,-0.2,-6,8,510,"");
  // axis1->SetName("axis1");
  // axis1->Draw();

  // TGaxis *axis2 = new TGaxis(-5,0.2,6,0.2,0.001,10000,510,"G");
  // axis2->SetName("axis2");
  // axis2->Draw();

//   TGaxis *axis3 = new TGaxis(-9,-0.8,-9,0.8,-8,8,50510,"");
//   axis3->SetName("axis3");
//   axis3->Draw();

//   TGaxis *axis4 = new TGaxis(-7,-0.8,-7,0.8,1,10000,50510,"G");
//   axis4->SetName("axis4");
//   axis4->Draw();

//   TGaxis *axis5 = new TGaxis(-5,-0.6,6,-0.6,1.2,1.32,80506,"-+");
//   axis5->SetName("axis5");
//   axis5->SetLabelSize(0.03);
//   axis5->SetTextFont(72);
//   axis5->SetLabelOffset(0.025);

//   axis5->Draw();

//   TGaxis *axis6 = new TGaxis(-5,0.6,6,0.6,100,900,50510,"-");
//   axis6->SetName("axis6");
//   axis6->Draw();

//   TGaxis *axis7 = new TGaxis(8,-0.8,8,0.8,0,9000,50510,"+L");
//   axis7->SetName("axis7");
//   axis7->SetLabelOffset(0.01);
//   axis7->Draw();

   TCanvas *c2 = new TCanvas("c2","c2",10,10,700,500);
   gStyle->SetOptStat(0);
   TH2F *h2 = new TH2F("h","Axes",100,0,10,100,-2,2);
   h2->Draw();
   TF1 *f1=new TF1("f1","-x",-10,10);
   TGaxis *A1 = new TGaxis(0,2,10,2,"f1",510,"-");
   A1->SetTitle("axis with decreasing values");
   A1->Draw();
   TF1 *f2=new TF1("f2","exp(x)",0,2);
   TGaxis *A2 = new TGaxis(1,1,9,1,"f2");
   A2->SetTitle("exponential axis");
   A2->SetLabelSize(0.03);
   A2->SetTitleSize(0.03);
   A2->SetTitleOffset(1.2);
   A2->Draw();
   TF1 *f3=new TF1("f3","log10(x)",1,1000);
   TGaxis *A3 = new TGaxis(2,-2,2,0,"f3",505,"G");
   A3->SetTitle("logarithmic axis");
   A3->SetLabelSize(0.03);
   A3->SetTitleSize(0.03);
   A3->SetTitleOffset(1.2);
   A3->Draw();
   return c2;
  
}


