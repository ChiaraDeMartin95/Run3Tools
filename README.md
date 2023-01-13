# Run3Tools

PostProcessing_Filters.C -> takes the output of the filter task (https://github.com/AliceO2Group/O2Physics/blob/master/EventFiltering/PWGLF/strangenessFilter.cxx) and post processes it 

Yileds_from_invmass.C -> takes mass vs pt histos and projects them along the mass axis in different pt intervals. Then, it fits invariant mass distributions with one or two gaussian + bkg function (pol1 or pol2) and produces summaty histos (mean, sigma, purity and yield per event vs pt)
