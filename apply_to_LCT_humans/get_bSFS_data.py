import sys, os
import numpy as np
import itertools
import scipy.special
import math
from tqdm import tqdm

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

def get_split_array_idxs(positions_of_muts, blocksize,num_blocks):
	intervals = np.arange(0, stop=num_blocks*blocksize+1, step=blocksize)
	bin_pos = np.digitize(positions_of_muts,intervals)
	bin_count = np.bincount(bin_pos)[1:]
	cs_bin_count = np.cumsum(bin_count)
	return cs_bin_count

def vectorize_all_blocks(sequence, num_mutypes):
	return [vectorize_all_combos(subseq, num_mutypes) for subseq in tqdm(sequence)]

def vectorize_all_combos(sequence, num_mutypes):
	if sequence.shape[-1]>0:
		sort_sequence = np.sort(sequence, axis=1)
		unique, counts = np.unique(sort_sequence, axis=0, return_counts=True)
		bSFS = np.apply_along_axis(to_bSFS, 1, unique, [num_mutypes,])
		return (bSFS, counts)
	else:
		return (np.zeros((1,num_mutypes), dtype=np.int8), np.array([sequence.shape[0],]))

def to_bSFS(row, args):
	num_mutypes, = args
	#result includes counts of zero at first
	#we can also encounter mutations ancestral to all samples
	result = np.zeros(num_mutypes+1, dtype=np.int8)
	unique, counts = np.unique(row, return_counts=True)
	result[unique] = counts
	#we are not interested in counting zeros
	return result[1:]

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

def comb_index(n, k):
    count = scipy.special.comb(n, k, exact=True)
    index = np.fromiter(itertools.chain.from_iterable(itertools.combinations(range(n), k)), 
                        int, count=count*k)
    return index.reshape(-1, k)

def generate_all_combinations(genotype_matrix, num_subsamples):
	num_lineages = genotype_matrix.shape[1]
	idx = comb_index(num_lineages, num_subsamples)
	g = genotype_matrix[:,idx]
	g = np.swapaxes(g, 0, 1)
	return g

def process_genotype_matrix(genotype_matrix, num_subsamples):
	#create new shape (comparison, position, genotype1, genotype2)
	all_combinations_genotype_matrix = generate_all_combinations(genotype_matrix, 2)
	print('made all combinations')
	#reshape allel genotype matrix 
	all_combinations_genotype_matrix = np.reshape(all_combinations_genotype_matrix, all_combinations_genotype_matrix.shape[:-2]+(-1,))	
	mapped_mutypes = map_mutypes(all_combinations_genotype_matrix, fully_labeled=False)
	print('mutypes mapped')
	return mapped_mutypes

def adapt_matrix_to_blocks(all_positions, pos_sweep, blocksize, num_blocks_tuple):
	num_blocks_left, num_blocks_right = num_blocks_tuple
	cutoff_left = np.argmin(all_positions<=pos_sweep - num_blocks_left * blocksize)
	cutoff_right = np.argmax(all_positions>pos_sweep + num_blocks_right * blocksize)
	return (cutoff_left, cutoff_right)

def get_num_blocks(all_positions, pos_sweep, blocksize):
	begin, end = all_positions[0], all_positions[-1]
	num_blocks_left = int(np.floor((pos_sweep - begin)/blocksize))
	num_blocks_right = int(np.floor((end - pos_sweep)/blocksize))
	return (num_blocks_left, num_blocks_right)

def mapped_mutypes_to_bsfs(mapped_mutypes, blocksize, num_blocks, num_subsamples, positions):
	split_array_idxs = get_split_array_idxs(positions, blocksize, num_blocks)
	mapped_mutypes_per_site = np.array_split(mapped_mutypes, split_array_idxs, axis=1)
	print('blocks made')
	to_bSFS = vectorize_all_blocks(mapped_mutypes_per_site[:-1], num_subsamples-1)
	print('vectorized all blocks')
	return to_bSFS

def main():
	directory = '' #containing filtered numpy array
	out_directory = directory
	blocksize = 1000
	name = f'lactase_all_individuals_{blocksize}' 
	num_subsamples = 4
	fully_labeled = False

	#pos_sweep = 136608646-134608673 
	pos_sweep =  1999973 #after setting left most to 0
	genotype_matrix = np.load(os.path.join(directory,'genotype_matrix_all_individuals.npy'))
	num_lineages = genotype_matrix.shape[1]
	positions = np.load(os.path.join(directory,'pos_all_individuals_norm.npy'))
	num_blocks_tuple = get_num_blocks(positions, pos_sweep, blocksize)
	num_blocks = sum(num_blocks_tuple)
	print(num_blocks)
	cutoff_left, cutoff_right = adapt_matrix_to_blocks(positions, pos_sweep, blocksize, num_blocks_tuple)
	genotype_matrix = genotype_matrix[cutoff_left:cutoff_right]
	positions = positions[cutoff_left:cutoff_right]
	mapped_mutypes = process_genotype_matrix(genotype_matrix, num_subsamples)
	result_list = mapped_mutypes_to_bsfs(mapped_mutypes, blocksize, num_blocks, num_subsamples, positions)
	process_results(result_list, out_directory, name, num_blocks, num_lineages, num_subsamples, num_subsamples-1)

if __name__ == "__main__":
	main() 
