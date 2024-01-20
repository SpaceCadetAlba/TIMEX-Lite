# TIMEX-Lite
TIMEX-Lite Onset Detection

This repo contains the MATLAB script for the TIMEX-Lite onset detection method.
This is a MATLAB script for musical onset detection associated with this paper: 10.17743/jaes.2021.0040
The method itself is derived based on the TIMEX onset detection method detailed in this paper: 10.1080/14015439.2018.1452977

The method itself is tags onsets by peakfinding from pitch slope.
Specifically this script finds onsets of the "Ta" sound, from the "T" noise burst.
The purpose of the onset detection method is for synchrony analysis of vocal musical performance tasks where each note uses the "Ta" sound.
By controlling performance task and onset detection in this way - we can achieve high precision synchrony analysis - this is useful for evaluating synchrony in contexts such as Network Music Performance.

Some tidying up of this repo needs done:
- Pack code blocks into functions.
- Add licence.

~PC
