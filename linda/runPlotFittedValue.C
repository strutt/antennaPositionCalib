void runPlotFittedValue() {
  // gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libAnitaEvent.so");
  //  gSystem->CompileMacro("plotFittedValueHPOL.C","k");
  gSystem->CompileMacro("plotFittedValueHPOL2.C","k");
  //  gSystem->CompileMacro("plotFittedValue.C","k");

  plotFittedValueHPOL2();
  
}
