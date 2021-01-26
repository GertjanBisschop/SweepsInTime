import msprime, tskit
import multiprocessing
import numpy as np
from tqdm import tqdm

def simulate_replicate(seed):
	sample_size = 24
	ts = msprime.simulate(sample_size, Ne=10000, recombination_rate=1e-7, mutation_rate=1.25e-7,random_seed=seed, length=1000000)
	ts.dump(f'/scratch/gbisshop/classic_sweeps/neutral_n24/{seed}.trees')

def main():
	replicates = 10000
	threads = 60
	seeds = np.random.randint(1, 2 ** 32, replicates)
	with multiprocessing.Pool(processes=threads) as pool:
		result = list(tqdm(pool.imap(simulate_replicate, seeds), desc="sim progress", total=len(seeds)))
	#simulate_replicate(seeds[0])

if __name__ == "__main__":
	main()
