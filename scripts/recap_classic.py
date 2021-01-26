#!/usr/bin/env python3
import multiprocessing
import msprime, pyslim
import os, glob, sys
import numpy as np
from tqdm import tqdm

def recap_mutate(tree):
	directory_out='/scratch/gbisshop/classic_sweeps/selection_recap'
	directory = '/scratch/gbisshop/classic_sweeps/selection_trees'
	done_file_list = [file for file in os.listdir(directory_out) if file.endswith('.recap')]
	new_name = tree.rstrip('.trees')+'.recap'
	if new_name not in done_file_list:
		ts = pyslim.load(os.path.join(directory, tree))
		ts = ts.recapitate(Ne=10000, recombination_rate=1e-7).simplify()
		ts_mutated = msprime.mutate(ts, rate=1.25e-7)
		ts_mutated.dump(os.path.join(directory_out,tree.rstrip('.trees')+'.recap'))
	return None

def main():
	directory = '/scratch/gbisshop/classic_sweeps/selection_trees'
	directory_out = '/scratch/gbisshop/classic_sweeps/selection_recap'
	threads = 60

	file_list = [file for file in os.listdir(directory) if file.endswith('.trees')]
	with multiprocessing.Pool(processes=threads) as pool:
		pool.map(recap_mutate, file_list)
	#for file in tqdm(file_list):
	#	recap_mutate(file)
print('done')

if __name__ == '__main__':
	main()
