{

  AnitaGeomTool* geom = AnitaGeomTool::Instance();

  Int_t ant1 = 0;
  Int_t ant2 = 16;

  // How degenerate a problem is this?
  // For a given theta/phi plane wave arrival, calculate deltaTExpected varying r+phi of antenna pair.



  // If you only had one event... with one pair of antennas, you would have 6 free parameters 2*(r,phi,z)
  // With only one deltaT.

  // Probably the deltaT picks out a 6d "contour" in the 6d parameter space.
  // Then with the next event at a slightly different source theta/phi you pick out a different contour.
  // You would then overlay all those contours and hopefully a small region would be favored by all those
  // fits.


  // Another way to pose the problem.
  // Suppose you mapped out the real deltaT as a function of theta and phi.
  // The question is what set of 2*(r,z,phi) best fits that parameter space?
  // You then run your favourite fitter with 6 free parameters.
  // Can you do that independently for each pair of antennas?
  // You can... but will the answers all converge?



}
