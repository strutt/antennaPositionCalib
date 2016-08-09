void runPlotAngResPos() {
  // gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  //  gSystem->AddIncludePath("-I${PLOTTER_DIR}");
  gSystem->AddIncludePath("-I/home/lindac/ANITA/Software/Utils/include");
  			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  //  gSystem->Load("/usr/lib/libfftw3.so");
  gSystem->Load("libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaCorrelator.so");
  gSystem->CompileMacro("plotAngResPos.C","k");

  // plotAngResPos("LDBVPOL", "LDBVPOL");
  //  plotAngResPos("WAIS", "LDBHPOL");
  plotAngResPos("WAIS", "WAIS");
  //  plotAngResPos("WAIS", "WAIS");

}
