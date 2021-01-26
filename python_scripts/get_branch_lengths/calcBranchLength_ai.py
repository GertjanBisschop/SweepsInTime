#!/usr/bin/env python3

import tskit
import os, glob, sys
import numpy as np

def main():

    directory = '' #set directory containing trees
    outputdir = ''


    file_list = [file for file in os.listdir(directory) if file.endswith('.recap')]

    Ne = 10000
    sample_size = 2
    #create numpy array
    sites = [25000+100*i for i in range(250)]
    #sites=[1316,3289,6578,9867,13156]
    #results = [replicate, site, SFS]

    #for the different times after sweep completion
    sTime = [0,5000,10000,15000,17500,19000,20000]
        
    for T in sTime: 
                sampleT = T
                results = np.zeros((len(file_list), len(sites), sample_size-1))

                for idx, treefile in enumerate(file_list):
                        #simplify tree
                        ts = tskit.load(os.path.join(directory, treefile))

                        sample = list(np.random.choice([u for u in ts.samples() if ts.get_time(u)==sampleT] , 
                                size=sample_size, replace=False))
                
                        tss = ts.simplify(samples=sample)
                        
                        temp = tss.allele_frequency_spectrum(
                                mode="branch", polarised=True, span_normalise=True, windows="trees")
                                
                        for idx2, site in enumerate(sites):
                                #what if we are doing things twice for the same tree??
                                temp_tree = tss.at(site)
                                tree_index = temp_tree.get_index()

                                results[idx, idx2] = temp[tree_index][1:-1]

                #numpy save array
                results.astype('float32').tofile(outputdir+'/np_s{}_{}.dat'.
                        format(str(sample_size), str(T).replace('.', '_')))     
        
if __name__ == '__main__':
        main()
