# Run3Tools

PostProcessing_Filters.C -> takes the output of the filter task (https://github.com/AliceO2Group/O2Physics/blob/master/EventFiltering/PWGLF/strangenessFilter.cxx) and post processes it 

Yields_from_invmass.C -> takes mass vs pt histos and projects them along the mass axis in different pt intervals. Then, it fits invariant mass distributions with one or two gaussian + bkg function (pol1 or pol2) and produces summary histos (mean, sigma, purity and yield per event vs pt). It can take in input whaterver root file contains the necessary histograms (for the time being, the output of the filter task setting isFilter to 1, or the output of the v0cascadesqa task, or the output of the postprocessing task if isPostProcess is chosen)

CompareSigmaWidthPurity.C --> takes in input mean or sigma or purity or raw yield vs pt and compares two different files (e.g. data vs MC or different reconstruction passes). For the time being, it takes in input the output of the post processing macro (all V0s and cascades) or the output of Yields_from_invmass.C (just one particle at a time). It also (only) computes the pseudoefficiency, if desired.
