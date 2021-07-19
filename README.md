# SweepsInTime
Repository containing the supplementary figures as well as all notebooks and scripts for Bisschop et al. 2020. Sweeps In Time.

# Abstract
Current methods of identifying positively selected regions of the genome are limited by their underlying model in two key ways: the model cannot account for the timing of the adaptive event and the analytic predictions are limited to single nucleotide polymorphisms.  
Here we develop a tractable method of describing the effect of positive selection on the genealogical histories in the surrounding genome, explicitly modeling both the timing and context of the adaptive event. In addition, our framework allows us to go beyond simple polymorphism data. We are able to leverage information contained in patterns of linked variants, and even with very small sample sizes, our analytic framework has high power to identify historically adaptive regions of the genome and to correctly infer both the time and strength of selection. Finally, we derived the marginal distribution of genealogical branch lengths at a locus affected by selection acting at a linked site. This provides a much-needed link between current theoretical models to recent advances in simulation procedures that have allowed researchers both to examine the evolution of genealogical histories at the level of full chromosomes and build methods that attempt to reconstruct full ancestries from genome sequence data.

## Supplementary_Figures
This pdf contains all supporting figures for the main text.

## SI_main
Mathematica notebook S1.  The working machinery to obtain the results is provided in the .nb file.  A .pdf version showing all function definitions and results is included. 

## python_scripts
This folder contains the scripts necessary to simulate the different scenarios, recapitulate the generated tree sequences and extract branch lengths/bSFS from those tree sequences.

## mathematica notebooks
Folder containing the machinery necessary to churn out composite likelihoods using the instantaneous sweep approximations and to estimate sweep parameters for simulated replicates.

## running the notebooks
Running 

	slim python_scripts/slimulations/classic_sweep.txt

will generate a bunch of `.trees` files containing treesequences that still need to be recapitated.

This can be done by running:

	python python_scripts/process_slimulations/recap.py

Next generate `.csv` files containing the bsfs for each of the `.recap` files containing the recapitated tree.

	python python_scripts/get_bSFS/make_bSFS_subsampling.py

The result of this workflow can be processed by the `likelihood_inference_run.nb` notebook.

dependencies: 

	slim, msprime, tskit, numpy, tqdm, pandas, scipy
	
