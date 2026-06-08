**Macros used for the QA of strangeness filters for the skimming of pp data collected by ALICE in Run 3 (updated at the end of Run 3 skimming):**
- Yields_from_invmass.C
- CompareSigmaWidthPurity.C
- QAplots_runbyrun.C
- QAplots_SkimmedUnskimmed.C

**Macros used to do the QA of weak decay vertices (V0s, casacdes) in new Monte Carlo productions or new data reconstructions in Run 3:**
- PostProcessV0AndCascQA_AO2D.C
- PostProcessV0AndCascQA_AO2D_New2026.C (for productions since May 2026)
These macros take in input the output of the QA task: https://github.com/AliceO2Group/O2Physics/blob/master/PWGLF/Tasks/QC/v0cascadesqa.cxx

# Description of main macros in this repo

PostProcessing_Filters.C -> takes the output of the filter task (https://github.com/AliceO2Group/O2Physics/blob/master/EventFiltering/PWGLF/strangenessFilter.cxx) and post processes it 

Yields_from_invmass.C -> takes mass vs pt histos and projects them along the mass axis in different pt intervals. Then, it fits invariant mass distributions with one or two gaussian + bkg function (pol1 or pol2) and produces summary histos (mean, sigma, purity and yield per event vs pt). It can take in input whaterver root file contains the necessary histograms (for the time being, the output of the filter task setting isFilter to 1, or the output of the v0cascadesqa task, or the output of the postprocessing task if isPostProcess is chosen)

CompareSigmaWidthPurity.C --> takes in input mean or sigma or purity or raw yield vs pt and compares two different files (e.g. data vs MC or different reconstruction passes). For the time being, it takes in input the output of the post processing macro (all V0s and cascades) or the output of Yields_from_invmass.C (just one particle at a time). It also (only) computes the pseudoefficiency, if desired.

QAplots_runbyrun.C --> compares the selectivity of the different filters, run by run

QAplots_SkimmedUnskimmed.C --> computes the ratio between the selectivity of filters in the skimmed sample vs the original unskimmed sample, run by run
