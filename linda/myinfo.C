{


  // Hack for ROOT 6...
  // Why is this necessary? Why is ROOT so silly?
  // See http://comments.gmane.org/gmane.comp.lang.c%2B%2B.root/16834
  if(strcmp(gStyle->GetName(), "benStyle")!=0){ 

    TStyle* benStyle = new TStyle("benStyle", "Ben's Beautiful Style");
  
    // Colors
  
    //set the background color to white
    //   benStyle->SetFillColor(10);
    //   benStyle->SetFrameFillColor(10);
    //   benStyle->SetCanvasColor(10);
    //   benStyle->SetPadColor(10);
    //   benStyle->SetTitleFillColor(0);
    //   benStyle->SetStatColor(10);
    benStyle->SetFillColor(10);
    benStyle->SetFrameFillColor(0);
    benStyle->SetCanvasColor(0);
    benStyle->SetPadColor(0);
    benStyle->SetTitleFillColor(0);
    benStyle->SetStatColor(0);
   
    //dont put a colored frame around the plots
    benStyle->SetFrameBorderMode(0);
    benStyle->SetCanvasBorderMode(0);
    benStyle->SetPadBorderMode(0);
   
    //use the primary color palette
   
    //  benStyle->SetPalette(1,0);
    //gStyle->Reset();
    //gStyle->SetPalette(1);
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor color;
    color.InitializeColors();
    color.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    // const int NRGBs = 3, NCont = 999;
    // gStyle->SetNumberContours(NCont);
    // Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
    // Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
    // Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
    // Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
    // TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont); 
   
    benStyle->SetNumberContours(NCont);

    //set the default line color for a histogram to be black
    benStyle->SetHistLineColor(kBlack);
    benStyle->SetHistLineWidth(3);
   
    //set the default line color for a fit function to be red
    benStyle->SetFuncColor(kRed);
   
    //make the axis labels black
    benStyle->SetLabelColor(kBlack,"xyz");
   
    //set the default title color to be black
    benStyle->SetTitleColor(kBlack);
   
    // Sizes
   
    //set the margins
    benStyle->SetPadBottomMargin(0.1);
    benStyle->SetPadTopMargin(0.1);
    benStyle->SetPadRightMargin(0.15);
    benStyle->SetPadLeftMargin(0.1);
   
    //set axis label and title text sizes
    benStyle->SetLabelSize(0.03,"xyz");
    benStyle->SetLabelOffset(0.02,"x");
    benStyle->SetTitleSize(0.04,"xyz");
    benStyle->SetTitleOffset(1,"x");
    benStyle->SetTitleOffset(1,"yz");
    benStyle->SetStatFontSize(0.03);
    benStyle->SetTextSize(0.03);
    benStyle->SetTitleBorderSize(0);
    benStyle->SetStatBorderSize(0);
    TGaxis::SetMaxDigits(3);
   
    //set line widths
    benStyle->SetHistLineWidth(1);
    benStyle->SetFrameLineWidth(1);
    benStyle->SetFuncWidth(1);

    //Set Legend boxes to be pretty
    //   benStyle->SetLegendFillColor(19);
    benStyle->SetLegendBorderSize(1);
   
    // Misc
   
    //align the titles to be centered
    //benStyle->SetTitleAlign(22);
   
    //set the number of divisions to show
    benStyle->SetNdivisions(506, "xy");
   
    //turn off xy grids
    benStyle->SetPadGridX(0);
    benStyle->SetPadGridY(0);
   
    //set the tick mark style
    benStyle->SetPadTickX(1);
    benStyle->SetPadTickY(1);
   
    //show the fit parameters in a box
    benStyle->SetOptFit(1111);
   
    //Turn off title
    //  benStyle->SetOptTitle(0); //RJN
   
    //turn off all other stats
    benStyle->SetOptStat(0000000);
   
    //marker settings
    benStyle->SetMarkerStyle(0);
    //benStyle->SetMarkerSize(0.9);
   
    // Fonts
   
    const int kBenFont = 42;
   
    benStyle->SetStatFont(kBenFont);
    benStyle->SetLabelFont(kBenFont,"xyz");
    benStyle->SetTitleFont(kBenFont,"xyz");
    benStyle->SetTextFont(kBenFont);
   
    //  done
    benStyle->cd();


    gROOT->ForceStyle();
    //   gStyle->ls();

    // load my macros?
    // gROOT->ProcessLine(".L sc.C");

    const char* anitaUtilInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
    TString incDir = TString::Format("%s/include", anitaUtilInstallDir);
    gInterpreter->AddIncludePath(incDir.Data());  

    gROOT->ProcessLine("#include <AnitaConventions.h>");  
    gROOT->ProcessLine("#include <RootTools.h>");
    gROOT->ProcessLine("#include <UsefulAdu5Pat.h>");  

    
  }

  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaCorrelator.so");  
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->Load("libMagicDisplay.so");      
  gSystem->Load("libBensAnitaTools.so");    
}
