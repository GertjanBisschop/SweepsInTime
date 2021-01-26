# SweepsInTime
Repository containing all notebooks and scripts for Bisschop et al. 2020. Sweeps In Time.

#Abstract
Current methods of identifying positively selected regions of the genome are limited by their underlying model in two key ways: the model cannot account for the timing of the adaptive event and the analytic predictions are limited to single nucleotide polymorphisms.  
Here we develop a tractable method of describing the effect of positive selection on the genealogical histories in the surrounding genome, explicitly modeling both the timing and context of the adaptive event. In addition, our framework allows us to go beyond simple polymorphism data. We are able to leverage information contained in patterns of linked variants, and even with very small sample sizes, our analytic framework has high power to identify historically adaptive regions of the genome and to correctly infer both the time and strength of selection. Finally, we derived the marginal distribution of genealogical branch lengths at a locus affected by selection acting at a linked site. This provides a much-needed link between current theoretical models to recent advances in simulation procedures that have allowed researchers both to examine the evolution of genealogical histories at the level of full chromosomes and build methods that attempt to reconstruct full ancestries from genome sequence data.

## SI_main
Mathematica notebooks S1 (classic sweeps) and S2 (adaptive introgression).  The working machinery to obtain the results is provided in the .nb files.  A .pdf version of each is provided that in which the function definitions are supressed to access the results more easily. *The simulations folder contains the simulated marginal branch length distribution results for the full range of parameters considered. 

