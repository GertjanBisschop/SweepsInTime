import sys, os
import numpy as np
import pandas as pd
import collections
import itertools
import tskit
from functools import partial
import multiprocessing
from tqdm import tqdm
import scipy.special
from timeit import default_timer as timer
import math

def comb_index(n, k):
    count = scipy.special.comb(n, k, exact=True)
    index = np.fromiter(itertools.chain.from_iterable(itertools.combinations(range(n), k)), 
                        int, count=count*k)
    return index.reshape(-1, k)

def generate_all_mutypes(num_lineages):
	binary_num =  np.array([''.join(str(el) for el in sublist) for sublist in itertools.product((0,1), repeat=num_lineages)])
	return {key:int(key, 2) for key in binary_num}

def map_mutypes(genotype_matrix,fully_labeled=True):
	num_lineages = genotype_matrix.shape[-1]
	if fully_labeled:
		bits = 2**np.arange(num_lineages)[::-1]
		max_bit = int(np.ones(num_lineages).dot(bits))
		return np.mod(genotype_matrix.dot(bits),max_bit)
	else:
		return np.mod(np.sum(genotype_matrix, axis=2),num_lineages)

def deal_with_finite_sites(pos_array):
    idxs = np.insert((np.diff(pos_array)==0).astype(np.bool), 0, False)
    if np.any(idxs): 
        # if there are duplicates, get new values by incrementing by one
        new_values = pos_array[idxs]+1
        # get non-duplicate values
        uniq_values = pos_array[~idxs]
        # insert new_values in non-duplicated values (required sorted)
        new_idxs = np.searchsorted(uniq_values, new_values)
        # recursive call
        return deal_with_finite_sites(np.sort(np.insert(uniq_values, new_idxs, new_values)))
    # if there are no duplicated values
    return pos_array

def process_genotype_matrix(ts, num_lineages, num_subsamples):
	sequence_length=ts.sequence_length
	positions_of_mutations = np.array([int(site.position) for site in ts.sites()])
	positions_of_mutations = deal_with_finite_sites(positions_of_mutations)
	genotype_matrix = ts.genotype_matrix()
	if any(p>=sequence_length for p in positions_of_mutations):
		#how much do we need to trim off:
		clip_at = np.argmax(positions_of_mutations>=sequence_length)
		positions_of_mutations=positions_of_mutations[:clip_at]
		#genotype matrix m x n where m is the number of sites
		genotype_matrix = genotype_matrix[:clip_at]
	idx = comb_index(num_lineages, num_subsamples)
	g = genotype_matrix[:,idx]
	g = np.swapaxes(g, 0, 1)
	return (g, positions_of_mutations)

def get_split_array_idxs(positions_of_muts, blocksize,num_blocks):
	intervals = np.arange(0, stop=num_blocks*blocksize+1, step=blocksize)
	bin_pos = np.digitize(positions_of_muts,intervals)
	bin_count = np.bincount(bin_pos)[1:]
	cs_bin_count = np.cumsum(bin_count)
	return cs_bin_count

def vectorize_all_blocks(sequence, num_mutypes):
	return [vectorize_all_combos(subseq, num_mutypes) for subseq in sequence]

def vectorize_all_combos(sequence, num_mutypes):
	if sequence.shape[-1]>0:
		sort_sequence = np.sort(sequence, axis=1)
		unique, counts = np.unique(sort_sequence, axis=0, return_counts=True)
		bSFS = np.apply_along_axis(to_bSFS, 1, unique, [num_mutypes,])
		return (bSFS, counts)
	else:
		return (np.zeros((1,num_mutypes-1), dtype=np.int8), np.array([sequence.shape[0],]))

def to_bSFS(row, args):
	num_mutypes, = args
	#result includes counts of zero at first
	#we can also encounter mutations ancestral to all samples
	result = np.zeros(num_mutypes, dtype=np.int8)
	unique, counts = np.unique(row, return_counts=True)
	result[unique] = counts
	#we are not interested in counting zeros
	return result[1:]

def read_ts_file(directory, file, num_lineages, time, num_subsamples):
	ts = tskit.load(os.path.join(directory, file))
	samples_time = [u for u in ts.samples() if ts.node(u).time == time]
	subsamples = np.random.choice(samples_time, size=num_lineages, replace=False) 
	#all_subsample_combos_iter = itertools.combinations(subsamples, num_subsamples)
	ts_simple = ts.simplify(subsamples)
	return ts_simple

def process_single_ts(name_and_file, directory, out_directory, time, num_lineages, num_subsamples, num_mutypes, blocksize, fully_labeled=True):
	name, file = name_and_file
	ts = read_ts_file(directory, file, num_lineages, time, num_subsamples)
	sequence_length=int(ts.sequence_length)
	num_blocks = int(sequence_length//blocksize)
	if sequence_length%blocksize:
		num_blocks = int(sequence_length//blocksize)
		sequence_length = int(num_blocks * blocksize)
		ts = ts.keep_intervals([[0, sequence_length]], simplify=False).rtrim()
	#process genotype matrix
	all_combinations_genotype_matrix, positions_of_muts = process_genotype_matrix(ts, num_lineages, num_subsamples)
	mapped_mutypes_per_site = map_mutypes(all_combinations_genotype_matrix, fully_labeled=fully_labeled)
	all_combinations_genotype_matrix = None
	split_array_idxs = get_split_array_idxs(positions_of_muts, blocksize,num_blocks)
	#split processed genotype matrix in blocks
	mapped_mutypes_per_site = np.array_split(mapped_mutypes_per_site, split_array_idxs, axis=1)
	to_bSFS = vectorize_all_blocks(mapped_mutypes_per_site[:-1], num_mutypes)
	process_results(to_bSFS, out_directory, name, num_blocks, num_lineages, num_subsamples, num_mutypes)
	return name

def process_results(result_list, out_directory, name, num_blocks, num_lineages, num_subsamples, num_mutypes):
	with open(os.path.join(out_directory, f'result_{name}.csv'), 'w') as file:
		for pos, (bSFS, counts) in enumerate(result_list):
			if pos<num_blocks:
				for obs, count in zip(bSFS,counts):
					mutype = ', '.join([str(o) for o in obs])
					print(f'{pos}, {mutype}, {count}', file=file)

		if pos<num_blocks-1:
			empty = np.zeros(num_mutypes-1, dtype=np.int32)
			empty_count = int(math.factorial(num_lineages)/(math.factorial(num_lineages-num_subsamples)*math.factorial(num_subsamples)))
			while pos<num_blocks-1:
				pos+=1
				mutype = ', '.join([str(o) for o in empty])
				print(f'{pos}, {mutype}, {empty_count}', file=file)
		
def main():
	directory = '' #dir containing trees
	out_directory = ''
	name = 'iton_subsampling'
	blocksize = 100
	time=30000
	cores = 60
	num_lineages = 20
	num_subsamples = 4
	fully_labeled = False

	file_list = [file for file in os.listdir(directory) if file.endswith('recap')]
	
	if fully_labeled:
		num_mutypes = 2**(num_subsamples)
	else:
		num_mutypes = num_subsamples + 1 #all zero and all one are included
	
	if cores>1:	
		process_single_specified = partial(
							process_single_ts, directory=directory, out_directory=out_directory, time=time, 
							num_lineages=num_lineages, num_subsamples=num_subsamples, num_mutypes=num_mutypes, blocksize=blocksize,
							fully_labeled=fully_labeled
							)
		with multiprocessing.Pool(processes=cores) as pool:
			result_list = list(tqdm(pool.imap(process_single_specified, enumerate(file_list)), total=len(file_list)))
	else:
		result_list = [process_single_ts((name, file), directory, out_directory, time, num_lineages, num_subsamples, num_mutypes, blocksize, fully_labeled) for name,file in tqdm(enumerate(file_list))]
		
if __name__ == "__main__":
	main() 
