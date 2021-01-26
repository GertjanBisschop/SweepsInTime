#!/usr/bin/env python3

import tskit
import os, glob, sys
import numpy as np
from multiprocessing import Pool



def processTime(T):
        sampleT = 2*2*Ne-2*Ne*T
        
        results = np.zeros((len(file_list), len(sites), sample_size-1))

        for idx, treefile in enumerate(file_list):
                #simplify tree
                ts = tskit.load(os.path.join(directory, treefile))

                #temp to find times used
                #tree = ts.first()
                #for u in range(0,100):
                #    print(u,tree.parent(u),tree.time(u),sep='\t')

                #break
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
        return(sampleT)                     


def main():
    global Ne, sample_size, sites, file_list, directory, outputdir
    directory = '/ceph/users/gbisshop/output/classic_sweeps/selection/selection_recap'
    outputdir = '/ceph/users/dsetter/gf_sel/classic_sweeps/marginals'


    file_list = [file for file in os.listdir(directory) if file.endswith('.recap')]
    #file_list = file_list[1:3]
    #file_list = file_list[1:500]
    Ne = 10000
    sample_size = 4
    #create numpy array
    sites = [500000+1000*i for i in range(250)]
    #sites=[1316,3289,6578,9867,13156]
    #results = [replicate, site, SFS]

    #for the different times after sweep completion
    #sTime = [0,5000,10000,15000,17500,19000,20000]
    #sTime = [0,1000,2500,5000,10000,15000]    
    sTime = [0, 0.1, 0.25, 0.5, 1, 1.5,2]
    #sTime = [int(2*2*Ne - 2*Ne*x) for x in sTime]
    print(sTime)
    print("processing")
    with Pool(len(sTime)) as p:
    		ret = p.map(processTime,sTime)
    		print(ret)

    
if __name__ == '__main__':
        main()
        print("done")
        quit()
