#############################################################################################################
# ANTENNA POSITION CALIB										    #
#############################################################################################################

Flow chart for how to find the antenna positions.

1.) With the most complete timing calibration do

    generateDeltaTTree

    This gives you a tree containing the time of maximum correlation for each antenna pair for the WAIS pulses.

2.) Turn off all pitch and roll in AnitaConventions and recompile all anita libs using the photogrammetry positions.

    generateAngularResolutionTree

    You will end up with quite a few angular resolution trees before we're through, so put this in folder called  photogrammetryPositionsWithNoPitchRoll or something like that.

3.) Now you're in a position to find the best pitch/roll/heading offsets. Lucky you.
    There is a degeneracy between tilts/offsets applied to the GPS position (ADU5) and the coordinate system of the antennas. (i.e. we could just apply a tilt to the antenna coordinates).
    But its a simple step that gets us closer to the answer.

4.) 