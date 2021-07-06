import pandas as pd
import numpy as np
import os, glob, sys
from tqdm import tqdm

def read_files(directory):
	all_files = os.listdir(directory)
	all_dfs = []
	for file in tqdm(all_files):
		df = pd.read_csv(os.path.join(directory, file), sep='\t')
		all_dfs.append(df)
	return all_dfs

def merge_files(all_dfs):
	return pd.concat(all_dfs)

def rescale_alpha(df):
	r = 1e-7
	Ne = 10000
	df['s'] = df['alpha'].apply(lambda x: r*np.log(2*Ne)/x)
	return df

def to_count_data(df, cut_bins, label_bins):
	groups = df.groupby(pd.cut(df.s, cut_bins))
	counts = np.array(groups.size())
	df = pd.DataFrame([(x,y) for x,y in zip(label_bins, np.array(counts))],columns=['s','count'])
	return df

def main():
	base_directory = ""
	directory = "n{}_{}Ne"
	times = [0.1, 0.5, 1.0, 2.0]
	sample_sizes = [4, 12, 20]
	bins_strong_selection = [0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2]
	new_bins = [(x+y)/2 for x,y in zip(bins_strong_selection[:-1], bins_strong_selection[1:])]
	bins_weak_selection = [x/10 for x in bins_strong_selection]
	for sample_size in sample_sizes:
		for time in times:
			file_directory = base_directory + directory.format(sample_size, time)
			result_df = merge_files(read_files(file_directory))
			result_df = rescale_alpha(result_df)
			result_df.to_csv(base_directory+directory.format(sample_size, time)+"_all.csv", index=False)
			binned_counts = to_count_data(result_df, new_bins, bins_strong_selection)
			binned_counts.to_csv(base_directory+directory.format(sample_size, time)+"_counts.csv", index=False)

if __name__ == "__main__":
	main()